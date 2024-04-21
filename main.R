# Main R Script to run everything


# Libraries ----
library(neon4cast)
library(dplyr)
library(lubridate)
library(ggplot2)
library(arrow)
library(glue)
library(readr)
library(tidyr)
library(tidyverse)
library(zoo)
library(rjags)
install.packages("daymetr")
library(daymetr)
library(purrr)
library(ncdf4)
devtools::install_github("EcoForecast/ecoforecastR",force=TRUE)
devtools::install_github('eco4cast/neon4cast')

# STEP 1A -- GET TAREGTS ----
##' Download Targets for Terrestrial 
##' @return data.frame in long format with days as rows, and time, site_id, variable, and observed as columns
download_targets <- function(){
  readr::read_csv("https://data.ecoforecast.org/neon4cast-targets/terrestrial_daily/terrestrial_daily-targets.csv.gz", guess_max = 1e6)
}

##' Download Site metadata
##' @return metadata dataframe
download_site_meta <- function(){
  site_data <- readr::read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv") 
  site_data %>% filter(as.integer(terrestrial) == 1)
}

# Run functions to collect data
target1 <- download_targets()       ## Y variables
sites <- unique(target1$site_id)
site_data  <- download_site_meta()

# STEP 1B -- GET NOAA ----
## Define sites ----
sites_sel <- c("BART", "OSBS", "KONZ", "SRER")

## Define dates ----
forecast_date <- Sys.Date()
noaa_date <- Sys.Date() - days(1)  #Need to use yesterday's NOAA forecast because today's is not available yet

# Define function ----
##' Download NOAA GEFS weather forecast
##' @param forecast_date start date of forecast
##' @return dataframe
download_met_forecast <- function(forecast_date){
  noaa_date <- forecast_date - lubridate::days(1)  #Need to use yesterday's NOAA forecast because today's is not available yet
  
  ## connect to data
  df_future <- neon4cast::noaa_stage2(start_date = as.character(noaa_date))
  
  ## filter available forecasts by date and variable
  met_future <- df_future |> 
    dplyr::filter(datetime >= lubridate::as_datetime(forecast_date), 
                  variable %in% c("air_pressure", "air_temperature",
                                  "relative_humidity", "precipitation_flux",
                                  "surface_downwelling_shortwave_flux_in_air",
                                  "surface_downwelling_longwave_flux_in_air")) |> 
    dplyr::collect()
  
  ## aggregate to daily
  #met_future <- met_future %>% 
  #  mutate(datetime = lubridate::as_date(datetime)) %>% 
  #  group_by(datetime, site_id, parameter, variable) |> 
  #  summarize(prediction = mean(prediction),.groups = "drop") |> 
  #    mutate(air_temperature = air_temperature - 273.15) |> 
  #  select(datetime, site_id, variable, prediction, parameter)
  
  #met_future <- met_future %>%
  # tidyr::pivot_longer(
  #  cols = 3, 
  # names_to="variable",
  #values_to = "prediction"
  #)
  
  return(met_future)
}

# Download future predictions of met ----
met_future <- download_met_forecast(forecast_date) 

# Select the sites we want 
met_future_sel <- met_future %>% filter(site_id %in% sites_sel)

##' append historical meteorological data into target file
df_past <- neon4cast::noaa_stage3() %>%
  dplyr::filter(site_id %in% sites_sel,
                lubridate::year(datetime) >= 2020,
                variable %in% c("air_pressure", "air_temperature",
                                "relative_humidity", "precipitation_flux",
                                "surface_downwelling_shortwave_flux_in_air",
                                "surface_downwelling_longwave_flux_in_air")) |> 
  dplyr::collect()

colnames(df_past)
colnames(met_future_sel)

merged_noaa <- rbind(met_future_sel, df_past) 

# Aggregate to daily ----
merged_noaa_daily <- merged_noaa %>%
  #mutate(datetime = lubridate::ymd_hms(datetime)) %>%
  mutate(
    year = year(datetime),
    month=month(datetime),
    day=day(datetime),
    prediction = ifelse(variable=="air_temperature", prediction-273.15, prediction)) %>%
  group_by(year, month, day, site_id, variable, parameter) %>%
  summarize(predict_daily = mean(prediction, na.rm=T), .groups = "drop") %>%
  mutate(date = ymd(paste0(year, "-", month, "-", day))) %>%
  select(-c(year, month, day))

# Step 2 - CALIBRATION ----


# Step 3 -- FORECAST ----
#Load params 
load("data_download_code/data/ensembleParameters.Rdata")

#Load forecast function
source("forecast_code/forecast_function.R") 

# Define forecast time period
start <- Sys.Date()-30
end <- Sys.Date()

time_points <- seq(start, end, "1 day")
site_ensemble <- vector("list", length(params))
site <- c("BART", "OSBS", "KONZ", "SRER")

#setup ensemble sizes and sample from the parameter spaces.
ens <- 1000

#grab met data.
met <- merged_noaa_daily %>% filter(
  date >= ymd(start) &
    date <= ymd(end)) %>%
  rename(pred_daily = predict_daily)

for (i in seq_along(params)) {
  #write parameters into ensembles
  ENS <- vector("list", ens)
  for (j in seq_along(ENS)) {
    #sample met data
    for (s in site){
      met_temp <- met %>% filter(site_id == s) %>% select(-site_id)
      ens.met <- met_temp[which(met_temp$parameter == sample(1:30, 1)),] 
      ENS[[j]] <- list(params = params[[i]]$params[j,],
                       met = ens.met,
                       x_ic = params[[i]]$predict[j])
    }}
  #forecast
  mu <- ENS %>% purrr::map(nee_forecast) %>% dplyr::bind_cols() %>% as.data.frame() %>% `colnames<-`(time_points)
  #store outputs
  site_ensemble[[i]] <- list(data = ENS, forecast = mu)
}



#time series plot.
#take BART for example
mu <- site_ensemble[[1]]$forecast
ci <- apply(mu,1,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale
plot(time_points, ci[3,], type="l", ylim=c(min(ci), max(ci)), xlab = "Date", ylab="NEE", main="NEE Forecasts -- BART")
ecoforecastR::ciEnvelope(time_points,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
lines(time_points, ci[1,])
lines(time_points, ci[2,], col=2)

# Plot OSBS
mu2 <- site_ensemble[[2]]$forecast
ci2 <- apply(mu2,1,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale
plot(time_points, ci2[3,], type="l", ylim=c(min(ci2), max(ci2)), xlab = "Date", ylab="NEE", main="NEE Forecasts -- OSBS")
ecoforecastR::ciEnvelope(time_points,ci2[1,],ci2[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
lines(time_points, ci2[1,])
lines(time_points, ci2[2,], col=2)

# Plot KONZ
mu3 <- site_ensemble[[3]]$forecast
ci3 <- apply(mu3,1,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale
plot(time_points, ci3[3,], type="l", ylim=c(min(ci3), max(ci3)), xlab = "Date", ylab="NEE", main="NEE Forecasts -- KONZ")
ecoforecastR::ciEnvelope(time_points,ci3[1,],ci3[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
lines(time_points, ci3[1,])
lines(time_points, ci3[2,], col=2)

# Plot SRER
mu4 <- site_ensemble[[4]]$forecast
ci4 <- apply(mu4,1,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale
plot(time_points, ci4[3,], type="l", ylim=c(min(ci4), max(ci4)), xlab = "Date", ylab="NEE", main="NEE Forecasts -- SRER")
ecoforecastR::ciEnvelope(time_points,ci4[1,],ci4[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
lines(time_points, ci4[1,])
lines(time_points, ci4[2,], col=2)


# STEP 4 -- KALMAN FILTER ----

##' Kalman Filter: Analysis step
##' @param  mu.f = Forecast mean (vector)
##' @param  P.f  = Forecast covariance (matrix)
##' @param  Y    = observations, with missing values as NAs) (vector)
##' @param  R    = observation error covariance (matrix)
##' @param  H    = observation matrix (maps observations to states)
KalmanAnalysis <- function(mu.f,P.f,Y,R,H,I){
  obs = !is.na(Y) ## which Y's were observed?
  if(any(obs)){
    H <- H[obs,]                                              ## observation matrix
    if (dim(t(H))[1] == 1) {
      K <- P.f %*% t(t(H)) %*% solve(H%*%P.f%*%t(t(H)) + R[obs,obs])  ## Kalman gain
    } else{
      K <- P.f %*% t(H) %*% solve(H%*%P.f%*%t(H) + R[obs,obs])  ## Kalman gain
    }
    mu.a <- mu.f + K%*%(Y[obs] - H %*% mu.f)                  ## update mean
    P.a <- (I - K %*% H)%*%P.f                                ## update covariance
    ## Note: Here's an alternative form that doesn't use the Kalman gain
    ## it is less efficient due to the larger number of matrix inversions (i.e. solve)
    ## P.a <- solve(t(H)%*%solve(R[obs,obs])%*%(H) + solve(P.f))                             
    ## mu.a <- P.a %*% (t(H)%*%solve(R[obs,obs])%*%Y[obs] + solve(P.f)%*%mu.f)
  } else {
    ##if there's no data, the posterior is the prior
    mu.a = mu.f
    P.a = P.f
  }
  return(list(mu.a=mu.a,P.a=P.a))
}

# Start data assimilation
# define period
nt <- length(time_points)
nee <- target1 %>% filter(variable=="nee")
nsite <- length(site)
site_ids <- site

# filter nee observations
nee.filter <- nee %>% dplyr::filter(datetime %in% time_points) 

# convert nee to a more generalized data frame
nee.mean <- nee.sd <- data.frame(matrix(NA, nrow = nsite, ncol = nt)) %>%
  `rownames<-`(site) %>%
  `colnames<-`(time_points)
## assume constant standard error for test run. ##
nee.sd[,] <- 0.2
for (id in seq_along(site)) {
  ind <- which(nee.filter$datetime %in% time_points & nee.filter$site_id == site[id])
  nee.mean[id, as.character(nee.filter$datetime[ind])] <- nee.filter$observation[ind]
}

# data assimilation
ens <- 1000
mu.f  = matrix(NA,nsite,nt+1)  ## forecast mean for time t
mu.a  = matrix(NA,nsite,nt)  ## analysis mean for time t
P.f  = array(NA,c(nsite,nsite,nt+1))  ## forecast variance for time t
P.a  = array(NA,c(nsite,nsite,nt))  ## analysis variance for time t
x_ic <- site_ensemble %>% 
  purrr::map('data') %>%
  unlist(recursive = F) %>%
  purrr::map('x_ic') %>%
  unlist %>%
  matrix(nrow = ens, ncol = nsite)
mu.f[,1] <- colMeans(x_ic)
P.f[,,1] <- cov(x_ic)
for (t in seq_along(time_points)) {
  # run forecast
  x_ic <- site_ensemble %>% 
    purrr::map('data') %>%
    unlist(recursive = F) %>%
    purrr::map2(t, nee_forecast) %>%
    unlist %>%
    matrix(nrow = ens, ncol = nsite) 
  
  # update mu.f and P.f
  mu.f[,t+1] <- colMeans(x_ic)
  P.f[,,t+1] <- cov(x_ic)
  Y <- nee.mean[,t]
  R <- diag(nee.sd[,t], nsite)
  H <- I <- diag(1, nsite)
  
  # run analysis
  analysis <- KalmanAnalysis(mu.f[,t+1], P.f[,,t+1], Y, R, H, I)
  mu.a[,t] <- analysis$mu.a
  P.a[,,t] <- analysis$P.a
  
  # update site_ensemble
  update <- mvtnorm::rmvnorm(ens, analysis$mu.a, analysis$P.a, method = "svd")
  
  # update initial conditions.
  for (i in seq_along(site_ids)) {
    for (j in 1:ens) {
      site_ensemble[[i]]$data[[j]]$x_ic <- update[j, i]
    }
  }
}

# plot
for(i in 1:nsite){
  ci = rbind(mu.a[i,]-1.96*sqrt(P.a[i,i,]),mu.a[i,]+1.96*sqrt(P.a[i,i,]))
  plot(time_points,mu.a[i,],ylim=range(ci,na.rm=TRUE),type='n',main=site_ids[i],xlab="Date",ylab="NEE")
  ecoforecastR::ciEnvelope(time_points,ci[1,],ci[2,],col="lightBlue")
  lines(time_points,mu.a[i,],col=4)
  points(time_points,nee.mean[i,],col="red")
  legend("topleft", lty=c(NA,1), pch=c("o", NA), col=c("red", "lightBlue"), legend = c("Observation", "Forecast"))
}

# STEP 5 -- SUBMIT FORECAST ----

##' Save forecast and metadata to file, submit forecast to EFI
##' @param forecast dataframe
##' @param team_info list, see example
##' @param submit boolean, should forecast be submitted to EFI challenge
submit_forecast <- function(forecast,team_info,submit=FALSE){
  
  #Forecast output file name in standards requires for Challenge.
  # csv.gz means that it will be compressed
  year <- year(Sys.Date())
  month <- month(Sys.Date())
  day <- day(Sys.Date())
  forecast_file <- paste0("terrestrial_daily-", year, "-", month,
                          "-", day, "-sustainseers.csv.gz")
  
  ## final format tweaks for submission
  forecast = forecast |> mutate(model_id = "sustainseers", family="ensemble") |>
    relocate(project_id,model_id,reference_datetime,datetime) |>
    relocate(parameter,.before = variable)|>
    relocate(family,.before = parameter)
  
  #Write csv to disk
  write_csv(forecast, forecast_file)
  
  #Confirm that output file meets standard for Challenge
  neon4cast::forecast_output_validator(forecast_file)
  
  # Generate metadata
  model_metadata = list(
    forecast = list(
      model_description = list(
        forecast_model_id =  "sustainseer", ## current git SHA
        name = "Dynamic linear model of net ecosystem exchange",
        type = "empirical",
        repository = "https://github.com/EcoForecast/SustainabilitySeers/"   ## put your REPO here *******************
      ),
      initial_conditions = list(
        status = "propogates"
      ),
      drivers = list(
        status = "propagates",
        complexity = 7, 
        propagation = list(
          type = "ensemble",
          size = 31)
      ),
      parameters = list(
        status = "data_driven",
        complexity = 2 # slope and intercept (per site)
      ),
      random_effects = list(
        status = "propagates"
      ),
      process_error = list(
        status = "propagates"
      ),
      obs_error = list(
        status = "propagates"
      )
    )
  )
  
  ## this function needs to be restored
  #metadata_file <- neon4cast::generate_metadata(forecast_file, team_info$team_list, model_metadata)
  
  if(submit){
    neon4cast::submit(forecast_file = forecast_file, ask = FALSE) #metadata = metadata_file,
  }
  
}

team_info <- list(
  team_list = list(
    list(
      team_name = "SustainabilitySeers",
      name = c("Dongchen Zhang", "Breanna van Loenen", "Tessa Keeney", "Katie Losada"),
      email = "bvanloen@bu.edu",
      institution = "Boston University"
    )
  ),
  SustainabilitySeers = "SustainabilitySeers"
)


#load("./data_download_code/data/site_ensemble.Rdata")

# get date for file name
library(lubridate)
year <- year(Sys.Date())
month <- month(Sys.Date())
day <- day(Sys.Date())

# reformat forecasts into submission form
final <- list("vector")
for (i in seq(site_ensemble)){
  forecast <- bind_rows(site_ensemble[[i]]$forecast)
  
  n <- ncol(forecast)-32
  
  start <- ymd(colnames(forecast)[31])+1
  end <- start+n
  times <- as.character(seq(start, end, "1 day"))
  colnames(forecast)[c(32:ncol(forecast))] <- c(times)
  
  forecastx <- forecast %>%
    pivot_longer(
      cols = everything(),
      names_to = "datetime",
      values_to = "prediction"
    ) %>%
    mutate(datetime = ymd(datetime)) %>%
    mutate(parameter = rep(c(1:1000), 31),
           site_id = i,
           variable="nee",
           duration = "P1D",
           project_id = "neon4cast",
           reference_datetime = min(datetime))
  
  final[[i]] <- forecastx
  
}

final_forecast <- bind_rows(final)
class(final_forecast$datetime)

final_forecast <- final_forecast %>% mutate(site_id = case_when(
  site_id==1 ~ "BART",
  site_id==2 ~ "OSBS",
  site_id==3 ~ "KONZ",
  site_id==4 ~ "SRER"
))

# save forecast
#wd <- "~forecast_code/output/"
#final_forecastcsv <- write.csv(final, paste0(wd, "terrestrial_daily-", year, "-", month,
#                                            "-", day, "-sustainseers.csv.gz")) # figure out how to make this into a format that can be submitted.

# Submit forecast
submit_forecast(final_forecast, team_info, submit = TRUE) # Assuming you want to submit the forecast immediately

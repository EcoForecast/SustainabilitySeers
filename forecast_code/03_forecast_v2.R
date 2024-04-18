library(purrr)
library(lubridate)
load("~/SustainabilitySeers/data_download_code/data/ensembleParameters.Rdata")
source("~/SustainabilitySeers/Data_Download_Functions/GEFS_download.R")

# Load covariate data ----
source("~/SustainabilitySeers/data_download_code/01_datatargetdownload.R") # NEE
source("~/SustainabilitySeers/data_download_code/01A_NOAA_datadownload.R") # weather

#Load forecast function
source("~/SustainabilitySeers/forecast_code/forecast_function.R") 

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

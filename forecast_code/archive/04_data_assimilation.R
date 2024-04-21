library(purrr)
load("data_download_code/data/ensembleParameters.Rdata")
source("Data_Download_Functions/GEFS_download.R")
#Get the data ----
nee <- readr::read_csv("https://sdsc.osn.xsede.org/bio230014-bucket01/challenges/targets/project_id=neon4cast/duration=PT30M/terrestrial_30min-targets.csv.gz", guess_max = 1e6)
site_ids <- params %>% purrr::map('site_id') %>% unlist
nsite <- length(site_ids)
nee <- nee %>% dplyr::filter(variable=="nee",
                      site_id %in% site_ids) %>%
  dplyr::mutate(year=lubridate::year(datetime),
         month=lubridate::month(datetime), 
         day=lubridate::day(datetime)) %>%
  dplyr::group_by(year, month, day, site_id) %>%
  dplyr::summarize(
    nee_daily = mean(observation, na.rm=T)
  ) %>%
  dplyr::mutate(date = lubridate::ymd(paste0(year, "-", month, "-", day)))

#function for the forecasts.
nee_forecast <- function(ensemble, t = NULL) {
  params <- ensemble$params
  #grab met
  met <- ensemble$met
  temp <- met$pred_daily[which(met$variable == "air_temperature")]
  precip <- met$pred_daily[which(met$variable == "precipitation_flux")]
  humid <- met$pred_daily[which(met$variable == "relative_humidity")]
  Pres <- met$pred_daily[which(met$variable == "air_pressure")]
  LW <- met$pred_daily[which(met$variable == "surface_downwelling_shortwave_flux_in_air")]
  SW <- met$pred_daily[which(met$variable == "surface_downwelling_longwave_flux_in_air")]
  #forecast
  x_ic <- ensemble$x_ic
  if (is.null(t)) {
    mu <- rep(NA, length(temp))
    mu[1] <- x_ic
    for (t in 2:length(temp)) {
      new_nee <- mu[t-1]  + 
        params["betaX"]*mu[t-1] + 
        params["betaIntercept"] + 
        params["betaTemp"]*temp[t] + 
        params["betaPrecip"]*precip[t] + 
        params["betahumid"]*humid[t] +
        params["betaSWFlux"]*SW[t] +
        params["betaPress"]*Pres[t] +
        params["betaLWFlux"]*LW[t]
      mu[t] <- rnorm(1, new_nee, 1/sqrt(params["tau_add"]))
    }
    return(mu)
  } else {
    new_nee <- x_ic  + 
      params["betaX"]*x_ic + 
      params["betaIntercept"] + 
      params["betaTemp"]*temp[t] + 
      params["betaPrecip"]*precip[t] + 
      params["betahumid"]*humid[t] +
      params["betaSWFlux"]*SW[t] +
      params["betaPress"]*Pres[t] +
      params["betaLWFlux"]*LW[t]
    return(rnorm(1, new_nee, 1/sqrt(params["tau_add"])))
  }
}
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
    K <- P.f %*% t(H) %*% solve(H%*%P.f%*%t(H) + R[obs,obs])  ## Kalman gain
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
time_points <- seq(as.Date("2024-01-01"), as.Date("2024-01-31"), "1 day")
nt <- length(time_points)
# filter nee observations
nee.filter <- nee %>% dplyr::filter(date %in% time_points)
# convert nee to a more generalized data frame
nee.mean <- nee.sd <- data.frame(matrix(NA, nrow = nsite, ncol = nt)) %>%
  `rownames<-`(site_ids) %>%
  `colnames<-`(time_points)
## assume constant standard error for test run. ##
nee.sd[,] <- 0.2
for (id in seq_along(site_ids)) {
  ind <- which(nee.filter$date %in% time_points & nee.filter$site_id == site_ids[id])
  nee.mean[id, as.character(nee.filter$date[ind])] <- nee.filter$nee_daily[ind]
}

# data assimilation
load("~/SustainabilitySeers/data_download_code/data/site_ensemble.Rdata")
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
  lines(time_points,nee.mean[i,])
  legend("topleft", lty=1, col=c("black", "lightBlue"), legend = c("Observation", "Analysis"))
}



library(purrr)
load("~/SustainabilitySeers/Data/calibration.Rdata")
source("~/SustainabilitySeers/Data_Download_Functions/GEFS_download.R")
#set initial value of NEE for the forecasts.
x_ic <- 0.5177017
#combine parameters from ef.out object.
params <- do.call(rbind, ef.out$params)
#we can also setup ensemble sizes and sample from the parameter spaces.
if (TRUE) {
  ens <- 100
  sample.inds <- sample(seq_along(params[,1]), ens)
  params <- params[sample.inds,]
}
#grab met data.
# met <- as.data.frame(ef.out$data$Xf)#we will need to sample from GEFs in the future.
# #remove na.
# na.inds <- which(is.na(met[,2]))
# met <- met[-na.inds,]

#accessing GEFS.
time_points <- seq(as.Date("2021-01-01"), as.Date("2021-03-31"), "1 day")
met_variables <- c("precipitation_flux", "air_temperature", "relative_humidity", "surface_downwelling_shortwave_flux_in_air")

met <- list()
print(paste0("Downloading GEFS weather forecasts from ", time_points[1], " to ", time_points[length(time_points)]))
pb <- utils::txtProgressBar(min = 0, max = length(time_points), style = 3)
for (i in seq_along(time_points)) {
  met[[i]] <- GEFS_download(date = time_points[i], site_name = "HARV", variables = met_variables, is.daily = T)
  utils::setTxtProgressBar(pb, i)
}
met <- do.call(rbind, met)

ENS <- vector("list", ens)
for (i in seq_along(ENS)) {
  #sample met data
  ens.met <- met[which(met$parameter == sample(1:31, 1)),]
  ENS[[i]] <- list(params = params[i,],
                   met = ens.met,
                   x_ic = x_ic)
}
#function for the forecasts.
nee_forecast <- function(ensemble) {
  params <- ensemble$params
  #grab met
  met <- ensemble$met
  temp <- met$pred_daily[which(met$variable == "air_temperature")]
  precip <- met$pred_daily[which(met$variable == "precipitation_flux")]
  humid <- met$pred_daily[which(met$variable == "relative_humidity")]
  #forecast
  x_ic <- ensemble$x_ic
  mu <- rep(NA, length(temp))
  mu[1] <- x_ic
  for (t in 2:length(temp)) {
    new_nee <- mu[t-1]  + 
      params["betaX"]*mu[t-1] + 
      params["betaIntercept"] + 
      params["betaTemp"]*temp[t] + 
      params["betaPrecip"]*precip[t] + 
      params["betahumid"]*humid[t]
    mu[t] <- new_nee
  }
  mu
}
#run forecasts.
mu <- ENS %>% purrr::map(nee_forecast) %>% dplyr::bind_cols() %>% as.data.frame() %>% `colnames<-`(time_points)
#time series plot.
ci <- apply(mu,1,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale
plot(time_points, ci[3,], type="l", ylim=c(min(ci), max(ci)), xlab = "Date", ylab="NEE", main="NEE Forecasts")
lines(time_points, ci[1,])
lines(time_points, ci[2,], col=2)

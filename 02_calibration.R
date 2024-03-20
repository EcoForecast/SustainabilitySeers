library(rjags)
library(daymetr)
library(purrr)
library(ncdf4)
devtools::install_github("EcoForecast/ecoforecastR",force=TRUE)
#define period.
time.points <- seq(as.Date("2021-01-01"), as.Date("2021-12-31"), "1 day")

#grab met drivers.
nc <- nc_open("~/SustainabilitySeers/Data/ERA5.2021.nc")
precipitation <- ncdf4::ncvar_get(nc, "precipitation_flux")
air_temp <- ncdf4::ncvar_get(nc, "air_temperature") - 273.1
sw <- ncdf4::ncvar_get(nc, "surface_downwelling_shortwave_flux_in_air")

#grab NEE
nee.df <- readr::read_csv("https://sdsc.osn.xsede.org/bio230014-bucket01/challenges/targets/project_id=neon4cast/duration=PT30M/terrestrial_30min-targets.csv.gz", guess_max = 1e6)
nee.df <- nee.df[which(nee.df$site_id == "HARV" & 
                         lubridate::year(nee.df$datetime) == "2021" & 
                         nee.df$variable == "nee"),]

#match time for nee
realtime <- as.POSIXct(nc$dim$time$vals*3600, origin = "2021-01-01", tz = "UTC")

#calculate daily average.
nee_daily <- temp_daily <- precip_daily <- sw_daily <- c()
days <- sort(unique(as.Date(realtime)))
for (i in seq_along(days)) {
  ind1 <- which(as.Date(nee.df$datetime) == days[i])
  ind2 <- which(as.Date(realtime) == days[i])
  nee_daily <- c(nee_daily, mean(nee.df$observation[ind1]))
  temp_daily <- c(temp_daily, mean(air_temp[ind2]))
  precip_daily <- c(precip_daily, mean(precipitation[ind2]))
  sw_daily <- c(sw_daily, mean(sw[ind2]))
}
#replace na
na.ind <- which(is.na(nee_daily))
nee_daily[na.ind] <- NA

## fit the model
data <- list(x_ic = mean(nee_daily, na.rm = T),
             tau_ic = 1/sd(nee_daily, na.rm = T),
             a_obs = 1,
             r_obs = 1,
             a_add = 1,
             r_add = 1,
             n = length(nee_daily),
             y = nee_daily,
             Precip = precip_daily,
             Temp = temp_daily,
             SW = sw_daily)

ef.out <- ecoforecastR::fit_dlm(model=list(obs="y",fixed="~ 1 + X + Temp + Precip + SW"), data)
out <- do.call(rbind, ef.out$predict)
x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(out[,x.cols],2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

#time series plot
plot(days, ci[3,], type="l", ylim = range(nee_daily, na.rm = T), col=1, xlab="Date", ylab="NEE")
lines(days, ci[1,], col=1)
points(days, nee_daily, col=2, pch=20)
legend("bottomleft", legend=c("CI", "NEE Observation"), lty=c(1,NA), col=c(1,2), pch=c(NA, 20))

#diagnostics
BGR <- gelman.plot(ef.out$params)

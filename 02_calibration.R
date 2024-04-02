library(rjags)
library(daymetr)
library(purrr)
library(ncdf4)
devtools::install_github("EcoForecast/ecoforecastR",force=TRUE)
devtools::install_github('eco4cast/neon4cast')

#Get the data ----
source("~/SustainabilitySeers/Data_Download_Functions/01_datatargetdownload.R")

gc()

nee.df <- target1 %>% filter(variable=="nee")
le.df <- target1 %>% filter(variable=="le")

#define period.
time.points <- seq(as.Date("2021-01-01"), as.Date("2021-12-31"), "1 day")
df_past <- df_past %>% 
  filter(lubridate::date(datetime) %in% time.points,
         site_id == "HARV")

met_all <- df_past#rbind(met_future, df_past)
met_all <- met_all %>% tidyr::pivot_wider(
  names_from=variable,
  values_from=prediction
)


#grab NEE for one site
nee.df.harv <- nee.df[which(nee.df$site_id == "HARV" & 
                         lubridate::year(nee.df$datetime) == "2021" & 
                         nee.df$variable == "nee"),]

#match time for nee
#realtime <- as.POSIXct(met_future$dim$time$vals*3600, origin = "2021-01-01", tz = "UTC")

#calculate daily average.
met_all <- met_all %>%
  mutate(year=year(datetime),
                     month=month(datetime), 
                     day=day(datetime)) %>%
  group_by(year, month, day, site_id) %>%
  summarize(
    temp_daily = mean(air_temperature, na.rm=T),
    precip_daily = mean(precipitation_flux, na.rm=T),
    humid_daily = mean(relative_humidity, na.rm=T),
    pressure_daily = mean(air_pressure, na.rm=T)
  ) %>%
  mutate(date = ymd(paste0(year, "-", month, "-", day)))

nee.df.daily <- nee.df.harv %>%
  mutate(year=year(datetime),
         month=month(datetime), 
         day=day(datetime)) %>%
  group_by(year, month, day, site_id) %>%
  summarize(
    nee_daily = mean(observation, na.rm=T)
  ) %>%
  mutate(date = ymd(paste0(year, "-", month, "-", day)))

merged_nee_met <- left_join(nee.df.daily, met_all, by=c("year", "month", "day", "site_id"),
                            relationship="many-to-one")
nee_daily <- precip_daily <- temp_daily <- humid_daily <- rep(NA, length(time.points))
non.na.inds <- which(time.points %in% nee.df.daily$date)

nee_daily[non.na.inds] <- merged_nee_met$nee_daily
precip_daily[non.na.inds] <- merged_nee_met$precip_daily
temp_daily[non.na.inds] <- merged_nee_met$temp_daily
humid_daily[non.na.inds] <- merged_nee_met$humid_daily

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
             humid = humid_daily)

ef.out <- ecoforecastR::fit_dlm(model=list(obs="y",fixed="~ 1 + X + Temp + Precip + humid"), data)
out <- do.call(rbind, ef.out$predict)
x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(out[,x.cols],2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

#time series plot
plot(time.points, ci[3,], type="l", ylim = range(nee_daily, na.rm = T), col=1, xlab="Date", ylab="NEE")
lines(time.points, ci[1,], col=1)
points(time.points, nee_daily, col=2, pch=20)
legend("bottomleft", legend=c("CI", "NEE Observation"), lty=c(1,NA), col=c(1,2), pch=c(NA, 20))

#diagnostics
BGR <- gelman.plot(ef.out$params)

plot(ef.out$params)

effectiveSize(ef.out)







# ------------------------------ #
# other code
nee_daily <- temp_daily <- precip_daily <- sw_daily <- c()
days <- sort(unique(as.Date(met_future$datetime)))
for (i in seq_along(days)) {
  ind1 <- which(as.Date(nee.df$datetime) == days[i])
  ind2 <- which(as.Date(met_future$datetime) == days[i])
  nee_daily <- c(nee_daily, mean(nee.df$observation[ind1]))
  temp_daily <- c(temp_daily, mean(air_temperature[ind2]))
  precip_daily <- c(precip_daily, mean(precipitation_flux[ind2]))
  sw_daily <- c(sw_daily, mean(sw[ind2]))
}

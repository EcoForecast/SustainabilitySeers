# Purpose: Plot NOAA data
source("data_download_code/01A_NOAA_datadownload.R")
library(ggplot2)

# Plot ---- 
# BART SDLF
merged_noaa_daily %>%
  filter(site_id=="BART") %>%
  filter(variable == "surface_downwelling_longwave_flux_in_air") %>%
  ggplot(aes(x=date, y=predict_daily, color=as.factor(parameter))) +
  geom_line() +
  ggtitle("Surface Downwelling Longwave Flux")

# BART RH
merged_noaa_daily %>%
  filter(site_id=="BART") %>%
  filter(variable == "relative_humidity") %>%
  ggplot(aes(x=date, y=predict_daily, color=as.factor(parameter))) +
  geom_line() +
  ggtitle("Relative Humidity")

# BART Temp 
merged_noaa_daily %>%
  filter(site_id=="BART") %>%
  filter(variable == "air_temperature") %>%
  ggplot(aes(x=date, y=predict_daily, color=as.factor(parameter))) +
  geom_line() +
  ggtitle("Air Temperature")

# Do all sites by variable ---
# Temp
merged_noaa_daily %>%
  filter(variable == "relative_humidity") %>%
  ggplot(aes(x=date, y=predict_daily, color=as.factor(parameter))) +
  geom_line() +
  ggtitle("Air Temperature") +
  facet_wrap(~site_id)

# RH
merged_noaa_daily %>%
  filter(variable == "air_temperature") %>%
  ggplot(aes(x=date, y=predict_daily, color=as.factor(parameter))) +
  geom_line() +
  ggtitle("Relative Humidity") +
  facet_wrap(~site_id)

# Pressure
merged_noaa_daily %>%
  filter(variable == "air_pressure") %>%
  ggplot(aes(x=date, y=predict_daily, color=as.factor(parameter))) +
  geom_line() +
  ggtitle("Air Pressure") +
  facet_wrap(~site_id) # SRER site is oddly very low


# Precip
merged_noaa_daily %>%
  filter(variable == "precipitation_flux") %>%
  ggplot(aes(x=date, y=predict_daily, color=as.factor(parameter))) +
  geom_line() +
  ggtitle("Precipitation") +
  facet_wrap(~site_id)

# SWF
merged_noaa_daily %>%
  filter(variable == "surface_downwelling_shortwave_flux_in_air") %>%
  ggplot(aes(x=date, y=predict_daily, color=as.factor(parameter))) +
  geom_line() +
  ggtitle("Shortwave Flux") +
  facet_wrap(~site_id)

# LWF
merged_noaa_daily %>%
  filter(variable == "surface_downwelling_longwave_flux_in_air") %>%
  ggplot(aes(x=date, y=predict_daily, color=as.factor(parameter))) +
  geom_line() +
  ggtitle("Longwave Flux") +
  facet_wrap(~site_id)




# ----------------- PREAMBLE ---------------------------- #
# Purpose: Pull NOAA data for terrestrial forecast
# Created: 2/28/2024
# Authors: Breanna van Loenen, Dongchen Zhang, Tessa Keeney, Katherine Losada
# Last updated by: ""
# ------------------------------------------------------- #

# 1. Load Libraries ----
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

# Define sites ----
sites_sel <- c("BART", "OSBS", "KONZ", "SRER")

# Define dates ----
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
df_past <- neon4cast::noaa_stage3(version = "v12", endpoint = "data.ecoforecast.org", verbose = TRUE) %>%
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

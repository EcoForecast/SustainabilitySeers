# Purpose: Grab NOAA data and make initial plots


# Get future NOAA data
forecast_date <- Sys.Date()
noaa_date <- Sys.Date() - days(1)  #Need to use yesterday's NOAA forecast because today's is not available yet

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
                                  "relative_humidity", "precipitation_flux")) |> 
    dplyr::collect()
  
  ## aggregate to daily
  met_future <- met_future %>% 
    mutate(datetime = lubridate::as_date(datetime)) %>% 
    group_by(datetime, site_id, parameter, variable) |> 
    summarize(prediction = mean(prediction),.groups = "drop") |> 
    #    mutate(air_temperature = air_temperature - 273.15) |> 
    select(datetime, site_id, variable, prediction, parameter)
  
  #met_future <- met_future %>%
  # tidyr::pivot_longer(
  #  cols = 3, 
  # names_to="variable",
  #values_to = "prediction"
  #)
  
  return(met_future)
}

#target     <- merge_met_past(target)   ## append met data (X) into target file
met_future <- download_met_forecast(forecast_date) 
unique(met_future$variable)

##' append historical meteorological data into target file
df_past <- neon4cast::noaa_stage3() %>%
  dplyr::filter(site_id %in% sites,
                lubridate::year(datetime) >= 2020,
                variable %in% c("air_pressure", "air_temperature",
                                "relative_humidity", "precipitation_flux")) |> 
  dplyr::collect()

#check for one site
d_site1 <- df_past |> 
  dplyr::filter(site_id == "BART") |>
  dplyr::collect()

df_past <- df_past %>%
  select(-c(family, reference_datetime))

merged_noaa <- rbind(met_future, df_past) 

## visual check of data
# Make time series of LE
target1 %>% 
  filter(variable=="le",
         site_id=="BART") %>%
  ggplot(aes(x = datetime, y = observation)) +
  geom_point() +
  #facet_grid(~site_id, scale ="free") +
  ggtitle("LE")

# Make time series of NEE
target1 %>% 
  filter(variable=="nee",
         site_id=="BART") %>%
  ggplot(aes(x = datetime, y = observation)) +
  geom_point() +
  #facet_grid(~site_id, scale ="free") +
  ggtitle("NEE")

#Make time series of all future NOAA
met_future %>%
  filter(site_id=="BART") %>%
  ggplot(aes(x=datetime, y=air_temperature)) +
  geom_point() +
  ggtitle("NOAA Air Temp, March - April 2024")

# Make timeseries of ALL past and future NOAA
unique(merged_noaa$variable)
merged_noaa %>%
  filter(site_id=="BART") %>%
  ggplot(aes(x=datetime, y=prediction)) +
  geom_point() +
  facet_wrap(~variable) +
  ggtitle("NOAA Weather Variables, 2020-2024")

# Make plots of NOAA covariate data
temp = d_site1[d_site1$variable == "air_temperature",]
press = d_site1[d_site1$variable == "air_pressure",]
rh  = d_site1[d_site1$variable == "relative_humidity",]
pcp  = d_site1[d_site1$variable == "preciptation_flux",]

plot(temp$datetime, temp$prediction, type='l',
     main = "BART Temperature (K)",
     xlab = "Date",
     ylab = "Temperature (K)")
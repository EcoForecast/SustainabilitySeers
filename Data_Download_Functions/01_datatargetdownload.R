library(neon4cast)
library(dplyr)
library(lubridate)
library(ggplot2)
library(arrow)
library(glue)
library(readr)


##' Download Targets for Terrestrial 
##' @return data.frame in long format with days as rows, and time, site_id, variable, and observed as columns
download_targets <- function(){
  readr::read_csv("https://sdsc.osn.xsede.org/bio230014-bucket01/challenges/targets/project_id=neon4cast/duration=PT30M/terrestrial_30min-targets.csv.gz", guess_max = 1e6)
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
                  variable == "air_temperature") |> 
    dplyr::collect()
  
  ## aggregate to daily
  met_future <- met_future %>% 
    mutate(datetime = lubridate::as_date(datetime)) %>% 
    group_by(datetime, site_id, parameter) |> 
    summarize(air_temperature = mean(prediction), .groups = "drop") |> 
    #    mutate(air_temperature = air_temperature - 273.15) |> 
    select(datetime, site_id, air_temperature, parameter)
  
  met_future <- met_future %>%
    tidyr::pivot_longer(
      cols = 3, 
      names_to="variable",
      values_to = "prediction"
    )
  
  return(met_future)
}

##' append historical meteorological data into target file
df_past <- neon4cast::noaa_stage3()

#check for one site
d_site1 <- df_past |> 
  dplyr::filter(site_id == "BART") |>
  dplyr::collect()

# Try for all sites
noaa_all <- df_past |> 
  dplyr::filter(site_id %in% sites) |>
  dplyr::collect() # this won't work on local hard drive, can't allocate 2,6Gb of data

noaa_all <- noaa_all %>%
  select(-c(family, reference_datetime)) 

merged_noaa <- rbind(met_future, noaa_all) 

#target     <- merge_met_past(target)   ## append met data (X) into target file
met_future <- download_met_forecast(forecast_date) 

## visual check of data
target1 %>% 
  filter(variable=="le",
         site_id=="BART") %>%
  ggplot(aes(x = datetime, y = observation)) +
  geom_point() +
  #facet_grid(~site_id, scale ="free") +
  ggtitle("LE")

target1 %>% 
  filter(variable=="nee",
         site_id=="BART") %>%
  ggplot(aes(x = datetime, y = observation)) +
  geom_point() +
  #facet_grid(~site_id, scale ="free") +
  ggtitle("NEE")

met_future %>%
  filter(site_id=="BART") %>%
  ggplot(aes(x=datetime, y=air_temperature)) +
  geom_point() +
  ggtitle("NOAA Air Temp")



# Make plots of NOAA covariate data
temp = d_site1[d_site1$variable == "air_temperature",]
press = d_site1[d_site1$variable == "air_pressure",]
rh  = d_site1[d_site1$variable == "relative_humidity",]
pcp  = d_site1[d_site1$variable == "preciptation_flux",]

plot(temp$datetime, temp$prediction, type='l',
     main = "BART Temperature (K)",
     xlab = "Date",
     ylab = "Temperature (K)")
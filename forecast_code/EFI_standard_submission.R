#devtools::install_github('eco4cast/neon4cast')

# Pull in forecast
source("~/SustainabilitySeers/forecast_code/04_data_assimilation_v2.R")

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
wd <- "~/SustainabilitySeers/forecast_code/output/"
final_forecastcsv <- write.csv(final, paste0(wd, "terrestrial_daily-", year, "-", month,
                                                  "-", day, "-sustainseers.csv.gz")) # figure out how to make this into a format that can be submitted.

# Submit forecast
submit_forecast(final_forecast, team_info, submit = TRUE) # Assuming you want to submit the forecast immediately

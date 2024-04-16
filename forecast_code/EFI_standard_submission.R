devtools::install_github('eco4cast/neon4cast')

time_points <- seq(as.Date("2024-01-01"), as.Date("2024-01-31"), "1 day")

##' Save forecast and metadata to file, submit forecast to EFI
##' @param forecast dataframe
##' @param team_info list, see example
##' @param submit boolean, should forecast be submitted to EFI challenge
submit_forecast <- function(forecast,team_info,submit=FALSE){

  #Forecast output file name in standards requires for Challenge.
  # csv.gz means that it will be compressed
  forecast_file <- paste0("terrestrial","-",min(time_points),"-",team_info$SustainabilitySeers,".csv.gz")

  ## final format tweaks for submission
  forecast = forecast |> mutate(model_id = team_info$SustainabilitySeers, family="ensemble") |>
    relocate(model_id,datetime) |>
    relocate(parameter,.before = variable) |>
    relocate(family,.before = parameter)

  #Write csv to disk
  write_csv(forecast, forecast_file)

  #Confirm that output file meets standard for Challenge
  neon4cast::forecast_output_validator(forecast_file)

  # Generate metadata
  model_metadata = list(
    forecast = list(
      model_description = list(
        forecast_model_id =  system("git rev-parse HEAD", intern=TRUE), ## current git SHA
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


load("~/SustainabilitySeers/Data/site_ensemble.Rdata")
site_ensemble # figure out how to make this into a format that can be submitted.


# Submit forecast
submit_forecast(combined_forecasts, team_info, submit = FALSE) # Assuming you want to submit the forecast immediately

devtools::install_github('eco4cast/neon4cast')

##' Save forecast and metadata to file, submit forecast to EFI
##' @param forecast dataframe
##' @param team_info list, see example
##' @param submit boolean, should forecast be submitted to EFI challenge
submit_forecast <- function(forecast,team_info,submit=FALSE){

  #Forecast output file name in standards requires for Challenge.
  # csv.gz means that it will be compressed
  forecast_file <- paste0("terrestrial","-",min(forecast$reference_datetime),"-",team_info$SustainabilitySeers,".csv.gz")

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
        name = "Air temperature to water temperature linear regression plus assume saturated oxygen",
        type = "empirical",
        repository = "https://github.com/ecoforecast/EF_Activities"   ## put your REPO here *******************
      ),
      initial_conditions = list(
        status = "absent"
      ),
      drivers = list(
        status = "propagates",
        complexity = 1, #Just air temperature
        propagation = list(
          type = "ensemble",
          size = 31)
      ),
      parameters = list(
        status = "data_driven",
        complexity = 2 # slope and intercept (per site)
      ),
      random_effects = list(
        status = "absent"
      ),
      process_error = list(
        status = "absent"
      ),
      obs_error = list(
        status = "absent"
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
      name = "Katie Losada",
      email = "your.email@example.com",
      institution = "Your Institution"
    )
  ),
  SustainabilitySeers = "SustainabilitySeers"
)


# Define time points and meteorological variables
time_points <- seq(as.Date("2024-01-01"), as.Date("2024-01-31"), "1 day")
met_variables <- c("precipitation_flux", "air_temperature", "air_pressure", "relative_humidity", "surface_downwelling_shortwave_flux_in_air", "surface_downwelling_longwave_flux_in_air")

# Define the sites
sites <- c("SRER", "OSBS", "BART", "KONZ")

# Initialize a list to store forecast data
all_forecasts <- list()

# Loop through each site
for (site in sites) {
  # Initialize a list to store data for this site
  met <- vector("list", length = length(time_points))

  print(paste0("Downloading GEFS weather forecasts from ",
               time_points[1], " to ",
               time_points[length(time_points)],
               " for ", site))

  # Progress bar
  pb <- utils::txtProgressBar(min = 0, max = length(time_points), style = 3)

  # Loop through each time point
  for (j in seq_along(time_points)) {
    met[[j]] <- GEFS_download(date = time_points[j], site_name = site, variables = met_variables, is.daily = TRUE)
    utils::setTxtProgressBar(pb, j)
  }

  # Append forecast to the list
  all_forecasts[[site]] <- met
}

# Combine all forecast data into a single dataframe
combined_forecasts <- bind_rows(all_forecasts)

# Submit forecast
submit_forecast(combined_forecasts, team_info, submit = FALSE) # Assuming you want to submit the forecast immediately

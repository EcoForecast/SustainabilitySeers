##' @title GEFS_download
##' @name  GEFS_download
##' @author Dongchen Zhang
##' 
##' @param date time when GEFS weather forecasts should be downloaded.
##' @param site_name the capitalized four-letter name for NEON sites.
##' @param variables the variable names to be extracted.
##' @details This function download GEFS weather forecasts of NEON.
##' 
##' @return It returns a data frame consisting of ensemble GEFS hourly weather forecasts.
##' @author Dongchen Zhang
##' ##' @examples
##' \dontrun{
##' dat <- GEFS_download("2024-03-06", "HARV", "air_temperature")
##' }
GEFS_download <- function(date, site_name, variables = NULL){
  #check packages.
  packages <- c("arrow", "purrr", "dplyr", "lubridate")
  packages.exist <- unlist(lapply(packages, require, character.only = TRUE))
  #install packages.
  for (i in seq_along(packages.exist)) {
    if (!packages.exist[i]) {
      install.packages(packages[i])
      require(packages[i])
    }
  }
  #set variables.
  if (is.null(variables)) {
    variables <- c("surface_downwelling_longwave_flux_in_air", "surface_downwelling_shortwave_flux_in_air","precipitation_flux",                       
                   "air_pressure", "relative_humidity", "air_temperature", "northward_wind", "eastward_wind")
  }
  bucket <- paste0("neon4cast-drivers/noaa/gefs-v12/stage2/parquet/0/", date)
  while ("try-error" %in% class(try(ds <- arrow::open_dataset(arrow::s3_bucket(bucket = bucket,
                                                                               endpoint_override = "data.ecoforecast.org",
                                                                               anonymous=TRUE))))) {
    date <- lubridate::as_date(date) - lubridate::days(1)#if data doesn't exist, try last day.
    bucket <- paste0("neon4cast-drivers/noaa/gefs-v12/stage2/parquet/0/", date)
  }
  df <- ds %>%
    filter(variable %in% variables,
           site_id == site_name) %>%
    collect() %>%
    mutate(datetime = lubridate::ymd_hms(datetime)) %>%
    filter(datetime >= lubridate::as_date(date))
  return(df)
}
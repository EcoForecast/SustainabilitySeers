##' @title NEON_download_CO2_concentration
##' @name  NEON_download_CO2_concentration
##' @author Dongchen Zhang
##' 
##' @param lat latitude.  
##' @param lon longitude.
##' @param start_date start date.
##' @param end_date end date.
##' @param store_dir path to where the downloaded files will be stored.
##' @details This function download and extract NEON CO2 concentration from NEON API.
##' 
##' @return It returns a data frame consisting of datetime and value.
##' @author Dongchen Zhang
##' ##' @examples
##' \dontrun{
##' dat <- NEON_download_CO2_concentration(lat = 40.3, lon = -105.5, start_date = "2020-01-01", end_date = "2020-12-31", store_dir = "~")
##' }
NEON_download_CO2_concentration <- function(lon, lat, start_date, end_date, store_dir){
  #check if folder exists.
  if (!file.exists(store_dir)) {
    dir.create(store_dir)
  }
  packages <- c("neonstore", "swfscMisc", "neonUtilities", "purrr", "dplyr")
  packages.exist <- unlist(lapply(packages, require, character.only = TRUE))
  #install packages.
  for (i in seq_along(packages.exist)) {
    if (!packages.exist[i]) {
      install.packages(packages[i])
      require(packages[i])
    }
  }
  #grab NEON sites info.
  neonsites <- neonstore::neon_sites(api = "https://data.neonscience.org/api/v0", .token = Sys.getenv("NEON_TOKEN"))
  neonsites <- dplyr::select(neonsites, "siteCode", "siteLatitude", "siteLongitude") #select for relevant columns
  #Find site name from coordinates using distance.
  dis <- swfscMisc::distance(lat1 = lat, lon1 = lon, lat2 = neonsites$siteLatitude, lon2 = neonsites$siteLongitude)
  sitename <- neonsites$siteCode[which.min(dis)]
  #download data.
  productID <- "DP4.00200.001"
  neonstore::neon_download(product = productID, 
                           dir = store_dir, 
                           table = NA, 
                           site = sitename, 
                           start_date = start_date, 
                           end_date = end_date, 
                           type = "basic",
                           api = "https://data.neonscience.org/api/v0")
  if ("try-error" %in% class(try(
    df <- neonstore::neon_read(table = "nsae-basic",
                               product = productID, 
                               start_date = start_date, 
                               end_date = end_date, 
                               dir = store_dir) %>% 
    # dplyr::mutate(collectDate = as.character(as.Date(collectDate))) %>% 
    dplyr::select(c("timeBgn", "data.fluxCo2.stor.flux"))
  ))) {
    return(NULL)
  } else {
    # df <- aggregate(.~timeBgn, df, mean, na.action = NULL, na.rm = TRUE)
    return(df)
  }
}
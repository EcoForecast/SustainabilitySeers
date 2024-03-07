##' @title MODIS_download
##' @name  MODIS_download
##' @author Dongchen Zhang
##' 
##' @param lat  latitude.  
##' @param lon longitude.
##' @param var variable names, can be "GPP", "LAI", or "LC".
##' @param start_date start date.
##' @param end_date end date.
##' @details This function download point data of GPP, LAI, and land cover from MODIS products.
##' 
##' @return It returns a list consisting of type of variable, value of interest, and date for each value.
##' @author Dongchen Zhang
##' ##' @examples
##' \dontrun{
##' dat <- MODIS_download(lat = 40.3, lon = -105.5, var = "GPP", start_date = "2020-01-01", end_date = "2020-12-31")
##' }

MODIS_download <- function(lat, lon, var, start_date, end_date){
  #check if the MODIS package exists.
  pkg.exist <- require(MODISTools)
  #if not, we will install it.
  if (!pkg.exist) {
    install.packages("MODISTools")
    require(MODISTools)
  }
  #grab product and band IDs for requested variable.
  if (var == "GPP") {
    product <- "MYD17A2HGF"
    band <- "Gpp_500m"
    scale <- 0.0001
  } else if (var == "LAI") {
    product <- "MCD15A3H"
    band <- "Lai_500m"
    scale <- 0.01
  } else if (var == "LC") {
    product <- "MCD12Q1"
    band <- "LC_Type1"
    LC.types <- list("1"="Evergreen Needleleaf Forests",
                     "2"="Evergreen Broadleaf Forests",
                     "3"="Deciduous Needleleaf Forests",
                     "4"="Deciduous Broadleaf Forests",
                     "5"="Mixed Forests",
                     "6"="Closed Shrublands",
                     "7"="Open Shrublands",
                     "8"="Woody Savannas",
                     "9"="Savannas",
                     "10"="Grasslands",
                     "11"="Permanent Wetlands",
                     "12"="Croplands",
                     "13"="Urban and Built-up Lands",
                     "14"="Cropland/Natural Vegetation Mosaics",
                     "15"="Permanent Snow and Ice",
                     "16"="Barren",
                     "17"="Water Bodies")
  }

  dat <- MODISTools::mt_subset(product = product,
                               lat = lat,
                               lon = lon,
                               band = band,
                               start = start_date,
                               end = end_date,
                               progress = FALSE)
  if (var == "LC") {
    return(list(var = var, data = unlist(LC.types[dat$value]), date = dat$calendar_date))
  } else {
    return(list(var = var, data = dat$value*scale, date = dat$calendar_date))
  }
}

# To retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  stop("Usage: Rscript script.R latitude longitude variable start_date end_date")
}

lat <- as.numeric(args[1])
lon <- as.numeric(args[2])
var <- args[3]
start_date <- args[4]
end_date <- args[5]

result <- MODIS_download(lat, lon, var, start_date, end_date)

output_file <- paste("MODIS_data_", Sys.Date(), ".csv", sep="")
write.table(result, file = output_file, row.names = FALSE, col.names=c("variable", "value", "date"))
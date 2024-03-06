# ----------------- PREAMBLE ---------------------------- #
# Purpose: Pull NEON/MODIS data for terrestrial forecast
# Created: 2/28/2024
# Author: Breanna van Loenen
# Last updated by: ""
# ------------------------------------------------------- #

# 1. Load Libraries ----
repos = "http://cran.us.r-project.org"
get.pkg <- function(pkg){
  loaded <- do.call("require",list(package=pkg))
  if(!loaded){
    print(paste("trying to install",pkg))
    install.packages(pkg,dependencies=TRUE,repos=repos)
    loaded <- do.call("require",list(package=pkg))
    if(loaded){
      print(paste(pkg,"installed and loaded"))
    } 
    else {
      stop(paste("could not install",pkg))
    }    
  }
}
get.pkg("RCurl")
get.pkg("readr")
get.pkg("XML")
get.pkg("arrow")
get.pkg("devtools")
get.pkg("MODISTools")
get.pkg("EML")
get.pkg("cronR")
get.pkg("tidyverse")

# 2. Pull Data ----
end_date <- Sys.Date()

## Variable 1: MODIS GPP ----
MODISTools::mt_products()
MODISTools::mt_bands(product="MOD17A2H") ## terrestrial GPP

# Load for the first time, starting 2020 (we may want to begin earlier)
#subset1 <- MODISTools::mt_subset(product = "MOD17A2H",
#                                band = "Gpp_500m",
#                                lat=42.5369,
#                                lon=-72.1727,
#                                start="2020-01-01",
#                                end=end_date,
#                                km_lr = 1,
#                                km_ab = 1,
#                               site_name = "HARV")
#save(subset1,file="GPP_file.RData")

# Now create the loop
GPP_file = "GPP_file.RData"
if(file.exists(GPP_file)){
  load(GPP_file)
  last_update <- max(GPP_file$end) 
} else {
  subset <- MODISTools::mt_subset(product = "MOD17A2H",
                                  band = "Gpp_500m",
                                  lat=42.5369,
                                  lon=-72.1727,
                                  start="2020-01-01",
                                  end=end_date,
                                  km_lr = 1,
                                  km_ab = 1,
                                  site_name = "HARV")
  save(subset,file=GPP_file)
}

head(subset)

# Grab updated data from last update to present day ----
updates <- MODISTools::mt_subset(product = "MOD17A2H",
                                 band = "Gpp_500m",
                                 start=last_update,
                                 end=end_date,
                                 lat=42.5369,
                                 lon=-72.1727,
                                 km_lr = 10,
                                 km_ab = 10, 
                                 site_name="HARV")

GPPall <- rbind(subset1, updates)
save(GPPall, file=GPP_file) 

## Variable 2: NEON Solar Radiation----
install.packages("neonUtilities")
library(neonUtilities)
library(lubridate)

#download data and all associated files initially to check structure
#PAR <- neonUtilities::loadByProduct(dpID = "DP1.00024.001", site="HARV", package="expanded", startdate="2020-01")
#head(PAR)
#PAR1m <- rbind(PAR$PARPAR_1min)
#save(PAR1m,file="photorad.RData")

#Make a loop 
if(file.exists("photorad.Rdata")){  
  load("photorad.Rdata")
  last_update1 <- max(lubridate::date(PAR1m$endDateTime))
} else{ 
  
  PAR <- neonUtilities::loadByProduct(dpID = "DP1.00024.001", site="HARV", package="expanded", startdate="2020-01")
  
  PAR1m <- rbind(PAR$PARPAR_1min)
  
  save(PAR1m,file="photorad.Rdata")
  
  #save ancillary files, since it includes a readme and error codes
  PAR_anc <- PAR[-3]
  save(PAR_anc, file="PAR_readme.RData")
  
} ## end download

# Now create grab updates
updates1 <- neonUtilities::loadByProduct(dpID = "DP1.00024.001", site="HARV", package="expanded", startdate=as.character(last_update1),
                                         include.provisional=TRUE)
updates1x <- rbind(updates1$PARPAR_1min) # updates of the 1min radiation data

allRadData <- rbind(PAR1m, updates1x)

save(updates1x, file=photorad.Rdata)


## Variable 3: NEON Precipitation----
#download data and all associated files initially to check structure
precipitation <- neonUtilities::loadByProduct(dpID = "DP1.00006.001", site="HARV", package="expanded", startdate="2020-01")
head(precipitation) # I think we only need primary precipitation for these purposes..
precip1 <- rbind(precipitation$PRIPRE_5min) # choosing the 5 min but we can change to 30 if this becomes too intensive
save(precip1,file="precip.RData")

#Make a loop 
if(file.exists("precip.RData")){  
  load("precip.RData")
  last_update2 <- max(lubridate::date(precip1$endDateTime))
} else{ 
  
  precipitation <- neonUtilities::loadByProduct(dpID = "DP1.00006.001", site="HARV", package="expanded", startdate="2020-01")
  
  precip1 <- rbind(precipitation$PRIPRE_5min) 
  
  save(precip1,file="precip.RData")
  
  #save ancillary files, since it includes a readme and error codes
  precip_anc <- precipitation[-4]
  save(precip_anc, file="precip_readme.RData")
  
} ## end download

# Now create grab updates
updates2 <- neonUtilities::loadByProduct(dpID = "DP1.00006.001", site="HARV", package="expanded", startdate=as.character(last_update2),
                                         include.provisional=TRUE)
updates2x <- rbind(updates2$PRIPRE_5min) # updates of the 1min radiation data

allPrecipData <- rbind(precip1, updates2x)

save(updates2x, file=precip.RData)



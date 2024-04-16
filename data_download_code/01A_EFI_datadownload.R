# ----------------- PREAMBLE ---------------------------- #
# Purpose: Pull NEON data for terrestrial forecast
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

# Define functions ----
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

# Collect data ----
# Run functions to collect data
target1 <- download_targets()       ## Y variables
sites <- unique(target1$site_id)
site_data  <- download_site_meta()

# Filter data ----
# Select only target sites BART, OSBS, KONZ, SRER
sites_sel <- c("BART", "OSBS", "KONZ", "SRER")
target_sub <- target1 %>% filter(site_id %in% sites_sel) 

# Now select only nee
target_sub_nee <- target_sub %>% filter(variable=="nee") %>%
  rename(nee = observation) %>%
  select(-variable)

# Plot data ----
## Attempt 1
ggplot(target_sub_nee, aes(x=ymd_hms(datetime), y=nee)) +
  geom_line() +
  xlab("Datetime") +
  ylab("NEE(umol CO2 m-2 s-1)") +
  ggtitle("Net Ecosystem Exchange, 30min") +
  facet_wrap(~site_id)

# data collection doesn't seem to start at SRER or OSBS until 2019, so let's begin our forecast later
target_sub_nee <- target_sub_nee %>% filter(year(datetime) >= 2020)


# Plot data 2 
ggplot(target_sub_nee, aes(x=ymd_hms(datetime), y=nee)) +
  geom_line() +
  xlab("Datetime (2020 - 2024)") +
  ylab("NEE(umol CO2 m-2 s-1)") +
  ggtitle("Net Ecosystem Exchange, 30min") +
  facet_wrap(~site_id)


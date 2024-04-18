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

source("~/SustainabilitySeers/data_download_code/01_datatargetdownload.R") # get data

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


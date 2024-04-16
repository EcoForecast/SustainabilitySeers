library(rjags)
library(daymetr)
library(purrr)
library(ncdf4)
devtools::install_github("EcoForecast/ecoforecastR",force=TRUE)
devtools::install_github('eco4cast/neon4cast')

#Get the data ----
source("~/SustainabilitySeers/Data/01A_EFI_datadownload.R")
source("~/SustainabilitySeers/Data/01B_NOAA_datadownload.R")


# Merge the target and covariate data
### Target data = target_sub_nee --> already filtered by site and time
### NOAA data = merged_noaa_daily --> already filtered by site and time

## First need daily aggregation of NEE
daily.nee <- target_sub_nee %>%
  mutate(year = year(datetime),
         month=month(datetime),
         day=day(datetime)) %>%
  group_by(year, month, day, site_id) %>%
  summarize(nee_daily = mean(nee, na.rm=T), .groups = "drop") %>%
  mutate(date = ymd(paste0(year, "-", month, "-", day))) %>%
  select(-c(year, month, day))


merged_noaa_daily2 <- merged_noaa_daily %>%
  pivot_wider(
    names_from = "variable",
    values_from= "predict_daily"
  )

merged_nee_met <- left_join(daily.nee, merged_noaa_daily2, by=c("date", "site_id"),
                            relationship="one-to-many")

# Can subset the time period for calibration -- just do a year (last year)
merged_nee_met_sub <- merged_nee_met %>% filter(lubridate::year(date)>= 2023)

# Investigate NA values for met
na <- merged_nee_met %>% filter(is.na(predict_daily))
na_dates <- na$date
na_met <- merged_noaa %>% filter(ymd(date) %in% na_dates) # none, so these must just be missing from NOAA?

# in future iterations, maybe we impute these values?

# Define variables ----
colnames(merged_nee_met_sub)[5:10] <- c("ap", "temp", "precip", "humid", "LWF", "SWF")

## fit the model
### BART ----
bart <- merged_nee_met_sub %>% filter(site_id=="BART")
nee_daily.bart <- bart$nee_daily
ap.bart <- bart$ap
temp.bart <- bart$temp
precip.bart <- bart$precip
humid.bart <- bart$humid
LWF.bart <- bart$LWF
SWF.bart <- bart$SWF

data.bart <- list(x_ic = mean(nee_daily.bart, na.rm = T),
             tau_ic = 1/sd(nee_daily.bart, na.rm = T),
             a_obs = 1,
             r_obs = 1,
             a_add = 1,
             r_add = 1,
             n = length(nee_daily.bart),
             y = nee_daily.bart,
             Press = ap.bart,
             Precip = precip.bart,
             Temp = temp.bart,
             humid = humid.bart,
             LWFlux = LWF.bart,
             SWFlux = SWF.bart)

# DLM BART ----
bart.out <- ecoforecastR::fit_dlm(model=list(obs="y",fixed="~ 1 + X + Temp + Precip + humid + Press + LWFlux + SWFlux"), 
                                  data.bart)
plot(bart.out$params)
BGR <- gelman.plot(bart.out)

bart.mat <- do.call(rbind, bart.out$predict)
effectiveSize(bart.mat)
cor(bart.mat)

burnin <- 500
bart.burn <- window(bart.out$params,start=burnin)  ## remove burn-in
plot(bart.burn)    

x.cols <- grep("^x",colnames(bart.mat)) ## grab all columns that start with the letter x
ci.bart <- apply(bart.mat[,x.cols],2,quantile,c(0.025,0.5,0.975))

time.bart <- bart$date

#time series plot
plot(time.bart, ci.bart[3,], type="l", ylim = range(nee_daily.bart, na.rm = T), col=1, xlab="Date", ylab="NEE")
lines(time.bart, ci.bart[1,], col=1)
points(time.bart, nee_daily.bart, col=2, pch=20)
legend("bottomleft", legend=c("CI", "NEE Observation"), lty=c(1,NA), col=c(1,2), pch=c(NA, 20))


### OSBS ----
osbs <- merged_nee_met_sub %>% filter(site_id=="OSBS")
nee_daily.osbs <- osbs$nee_daily
ap.osbs <- osbs$ap
temp.osbs <- osbs$temp
precip.osbs <- osbs$precip
humid.osbs <- osbs$humid
LWF.osbs <- osbs$LWF
SWF.osbs <- osbs$SWF

data.osbs <- list(x_ic = mean(nee_daily.osbs, na.rm = T),
                  tau_ic = 1/sd(nee_daily.osbs, na.rm = T),
                  a_obs = 1,
                  r_obs = 1,
                  a_add = 1,
                  r_add = 1,
                  n = length(nee_daily.osbs),
                  y = nee_daily.osbs,
                  Press = ap.osbs,
                  Precip = precip.osbs,
                  Temp = temp.osbs,
                  humid = humid.osbs,
                  LWFlux = LWF.osbs,
                  SWFlux = SWF.osbs)

# DLM OSBS ----
osbs.out <- ecoforecastR::fit_dlm(model=list(obs="y",fixed="~ 1 + X + Temp + Precip + humid + Press + LWFlux + SWFlux"), 
                                  data.osbs)
plot(osbs.out$params)
BGR <- gelman.plot(osbs.out)

#burnin <- 500
#osbs.burn <- window(osbs.out$params,start=burnin)  ## remove burn-in
#plot(osbs.burn) 

osbs.mat <- do.call(rbind, osbs.out$predict)
effectiveSize(osbs.mat)
cor(osbs.mat)

x.cols <- grep("^x",colnames(osbs.mat)) ## grab all columns that start with the letter x
ci.osbs <- apply(osbs.mat[,x.cols],2,quantile,c(0.025,0.5,0.975))

time.osbs <- osbs$date

#time series plot
plot(time.osbs, ci.osbs[3,], type="l", ylim = range(nee_daily.osbs, na.rm = T), col=1, xlab="Date", ylab="NEE")
lines(time.osbs, ci.osbs[1,], col=1)
points(time.osbs, nee_daily.osbs, col=2, pch=20)
legend("bottomleft", legend=c("CI", "NEE Observation"), lty=c(1,NA), col=c(1,2), pch=c(NA, 20))


### KONZ ----
konz <- merged_nee_met_sub %>% filter(site_id=="KONZ")
nee_daily.konz <- konz$nee_daily
ap.konz <- konz$ap
temp.konz <- konz$temp
precip.konz <- konz$precip
humid.konz <- konz$humid
LWF.konz <- konz$LWF
SWF.konz <- konz$SWF

data.konz <- list(x_ic = mean(nee_daily.konz, na.rm = T),
                  tau_ic = 1/sd(nee_daily.konz, na.rm = T),
                  a_obs = 1,
                  r_obs = 1,
                  a_add = 1,
                  r_add = 1,
                  n = length(nee_daily.konz),
                  y = nee_daily.konz,
                  Press = ap.konz,
                  Precip = precip.konz,
                  Temp = temp.konz,
                  humid = humid.konz,
                  LWFlux = LWF.konz,
                  SWFlux = SWF.konz)
# DLM KONZ ----
konz.out <- ecoforecastR::fit_dlm(model=list(obs="y",fixed="~ 1 + X + Temp + Precip + humid + Press + LWFlux + SWFlux"), 
                                  data.konz)
plot(konz.out$params)
BGR <- gelman.plot(konz.out)

konz.mat <- do.call(rbind, konz.out$predict)
effectiveSize(konz.mat)
cor(konz.mat)

x.cols <- grep("^x",colnames(konz.mat)) ## grab all columns that start with the letter x
ci.konz <- apply(konz.mat[,x.cols],2,quantile,c(0.025,0.5,0.975))

time.konz <- konz$date

#time series plot
plot(time.konz, ci.konz[3,], type="l", ylim = range(nee_daily.konz, na.rm = T), col=1, xlab="Date", ylab="NEE")
lines(time.konz, ci.konz[1,], col=1)
points(time.konz, nee_daily.konz, col=2, pch=20)
legend("bottomleft", legend=c("CI", "NEE Observation"), lty=c(1,NA), col=c(1,2), pch=c(NA, 20))


### SRER ----
srer <- merged_nee_met_sub %>% filter(site_id=="SRER")
nee_daily.srer <- srer$nee_daily
ap.srer <- srer$ap
temp.srer <- srer$temp
precip.srer <- srer$precip
humid.srer <- srer$humid
LWF.srer <- srer$LWF
SWF.srer <- srer$SWF

data.srer <- list(x_ic = mean(nee_daily.srer, na.rm = T),
                  tau_ic = 1/sd(nee_daily.srer, na.rm = T),
                  a_obs = 1,
                  r_obs = 1,
                  a_add = 1,
                  r_add = 1,
                  n = length(nee_daily.srer),
                  y = nee_daily.srer,
                  Press = ap.srer,
                  Precip = precip.srer,
                  Temp = temp.srer,
                  humid = humid.srer,
                  LWFlux = LWF.srer,
                  SWFlux = SWF.srer)


# DLM SRER ----
srer.out <- ecoforecastR::fit_dlm(model=list(obs="y",fixed="~ 1 + X + Temp + Precip + humid + Press + LWFlux + SWFlux"), 
                                  data.srer)
plot(srer.out$params)
BGR <- gelman.plot(srer.out)

srer.mat <- do.call(rbind, srer.out$predict)
effectiveSize(srer.mat)
cor(srer.mat)

x.cols <- grep("^x",colnames(srer.mat)) ## grab all columns that start with the letter x
ci.srer <- apply(srer.mat[,x.cols],2,quantile,c(0.025,0.5,0.975))

time.srer <- srer$date

#time series plot
plot(time.srer, ci.srer[3,], type="l", ylim = range(nee_daily.srer, na.rm = T), col=1, xlab="Date", ylab="NEE")
lines(time.srer, ci.srer[1,], col=1)
points(time.srer, nee_daily.srer, col=2, pch=20)
legend("bottomleft", legend=c("CI", "NEE Observation"), lty=c(1,NA), col=c(1,2), pch=c(NA, 20))

# Final outputs to use (all 5000 in length, can't sample() MCMC class, and sampling the matrix requires independent
# samples of predictions and parameters, which doesn't seem legit)
## BART : bart.out (params, predic)
## OSBS : osbs.out (params, predic)
## KONZ : konz.out (params, predic)
## SRER : srer.out (params, predic)





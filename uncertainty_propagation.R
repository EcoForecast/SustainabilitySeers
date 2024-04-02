library(rjags)
library(daymetr)
library(purrr)
library(ncdf4)
devtools::install_github("EcoForecast/ecoforecastR",force=TRUE)

#Get the data ----
source("~/SustainabilitySeers/Data_Download_Functions/01_datatargetdownload.R")

gc()

nee.df <- target1 %>% filter(variable=="nee")
le.df <- target1 %>% filter(variable=="le")

#define period.
time.points <- seq(as.Date("2021-01-01"), as.Date("2021-12-31"), "1 day")
df_past <- df_past %>% 
  filter(lubridate::date(datetime) %in% time.points)

met_all <- rbind(met_future, df_past)
met_all <- met_all %>% tidyr::pivot_wider(
  names_from=variable,
  values_from=prediction
)


#grab NEE for one site
nee.df.harv <- nee.df[which(nee.df$site_id == "HARV" & 
                              lubridate::year(nee.df$datetime) == "2021" & 
                              nee.df$variable == "nee"),]

#match time for nee
#realtime <- as.POSIXct(met_future$dim$time$vals*3600, origin = "2021-01-01", tz = "UTC")

#calculate daily average.
met_all <- met_all %>%
  mutate(year=year(datetime),
         month=month(datetime), 
         day=day(datetime)) %>%
  group_by(year, month, day, site_id) %>%
  summarize(
    temp_daily = mean(air_temperature, na.rm=T),
    precip_daily = mean(precipitation_flux, na.rm=T),
    humid_daily = mean(relative_humidity, na.rm=T),
    pressure_daily = mean(air_pressure, na.rm=T)
  ) %>%
  mutate(date = ymd(paste0(year, "-", month, "-", day)))

nee.df.daily <- nee.df %>%
  mutate(year=year(datetime),
         month=month(datetime), 
         day=day(datetime)) %>%
  group_by(year, month, day, site_id) %>%
  summarize(
    nee_daily = mean(observation, na.rm=T)
  ) %>%
  mutate(date = ymd(paste0(year, "-", month, "-", day)))

merged_nee_met <- left_join(nee.df, met_all, by=c("year", "month", "day", "site_id"),
                            relationship="many-to-one")

nee_daily <- merged_nee_met$nee_daily
precip_daily <- merged_nee_met$precip_daily
temp_daily <- merged_nee_met$temp_daily
humid_daily <- merged_nee_met$humid_daily

#replace na
na.ind <- which(is.na(nee.df.daily))
nee.df.daily[na.ind] <- NA

## fit the model
data <- list(x_ic = mean(nee_daily, na.rm = T),
             tau_ic = 1/sd(nee_daily, na.rm = T),
             a_obs = 1,
             r_obs = 1,
             a_add = 1,
             r_add = 1,
             n = length(nee_daily),
             y = nee_daily,
             Precip = precip_daily,
             Temp = temp_daily,
             humid = humid_daily)

ef.out <- ecoforecastR::fit_dlm(model=list(obs="y",fixed="~ 1 + X + Temp + Precip + humid"), data)
out <- do.call(rbind, ef.out$predict)
x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(out[,x.cols],2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

#time series plot
plot(days, ci[3,], type="l", ylim = range(nee_daily, na.rm = T), col=1, xlab="Date", ylab="NEE")
lines(days, ci[1,], col=1)
points(days, nee_daily, col=2, pch=20)
legend("bottomleft", legend=c("CI", "NEE Observation"), lty=c(1,NA), col=c(1,2), pch=c(NA, 20))

#diagnostics
BGR <- gelman.plot(ef.out$params)

plot(ef.out$params)

effectiveSize(ef.out)



#Forward Simulation

##settings
s <- 1
Nmc = 1000  #number of Monte Carlow draws
ylim = range(nee_daily, na.rm = T)   #not sure what to set limits as
N.cols <- c("black","red","green","blue","orange") ## set colors
trans <- 0.8       ## set transparancy
time = 1:(days*2)    ## total time
time1 = 1:days       ## calibration period
time2 = time1+days   ## forecast period

#select one site
plot.run <- function(){
  sel = seq(s,ncol(ci),by=31)
  plot(time,time,type='n',ylim=ylim,ylab="Y")
  ecoforecastR::fit_dlm(model=list(obs="Y",fixed="~ 1 + X + Temp + Precip"), data)
  lines(time1,ci[2,sel],col="blue")
  points(time1,No[s,])
}
  
ci <- apply(out[,x.cols],2,quantile,c(0.025,0.5,0.975))
plot.run()


forecastY <- function(x_ic,Precip,Temp,a_add=0,n=Nmc){
  Y <- matrix(NA,n,days)  ## storage
  Yprev <- x_ic          ## initialize
  for(t in 1:days){
    mu =    ## calculate mean (unsure what it should be)
    Y[,t] <- rlnorm(n,mu,Q)   ## predict next step
    Yprev <- Y[,t]        ## update IC
  }
  return(Y)
}



### Deterministic Prediction

## calculate mean of all inputs
precip.mean <- matrix(apply(precip_ensemble,2,mean),1,days) ## driver
temp.mean <- matrix(apply(temp_ensemble,2,mean),1,days)
## parameters
params <- as.matrix(out$params)
param.mean <- apply(params,2,mean)
## initial conditions
IC <- as.matrix(out$predict)

Y.det <- forecastY(x_ic=mean(nee_daily, na.rm = T),  #unsure of length 
                   precip=precip.mean,
                   temp = temp.mean,
                   beta=param.mean["beta"],
                   alpha=param.mean["alpha_site[6]"],     ##unsure what paramters should be
                   a_add=0,  ## process error off
                   n=1)

## Plot run
plot.run()
lines(time2,Y.det,col="purple",lwd=3)



##Monte Carlo Propogation

##sample parameter rows from previous analysis
prow = sample.int(nrow(params),Nmc,replace=TRUE)

Y.I <-  forecastY(IC=IC[prow,"Y[1,...]"],  ## sample IC (posterier distribution, site 6 and unsure what year to pick) 
                  precip=precip.mean,
                  temp=temp.mean,
                  beta=param.mean["beta"],  #unsure what parameters should be 
                  alpha=param.mean["alpha_site[6]"],
                  a_add=0,  #process error off
                  n=Nmc)

  #plot run 
plot.run()
Y.I.ci = apply(Y.I,2,quantile,c(0.025,0.5,0.975))
ecoforecastR::ciEnvelope(time2,Y.I.ci[1,],Y.I.ci[3,],col=col.alpha(N.cols[1],trans))
lines(time2,Y.I.ci[2,],lwd=0.5)



##Parameter Uncertainty 

Y.IP <- forecastY(IC=IC[prow,"N[6,...]"],  ## sample IC
                  precip=precip.mean,
                  temp = temp.mean, 
                  beta=params[prow,"beta"],     #unsure about parameters 
                  alpha=params[prow,"alpha_site[6]"],
                  a_add=0,
                  n=Nmc)

## plot run
plot.run()
Y.IP.ci = apply(Y.IP,2,quantile,c(0.025,0.5,0.975))
ecoforecastR::ciEnvelope(time2,Y.IP.ci[1,],Y.IP.ci[3,],col=col.alpha(N.cols[2],trans))
ecoforecastR::ciEnvelope(time2,Y.I.ci[1,],Y.I.ci[3,],col=col.alpha(N.cols[1],trans))
lines(time2,Y.I.ci[2,],lwd=0.5)




## Driver Uncertainty 


#forecast precipitation 
NE = 31    #NE is member ensemble of precipitation forecast 
#{r,echo=FALSE}
plot(time2,precip_ensemble[1,],type='n',ylim=range(precip_ensemble),xlab="time",ylab="precipitation (mm)")
for(i in 1:NE){  
  lines(time2,precip_ensemble[i,],lwd=0.5,col="grey")
}


#forecast temperature 
#{r,echo=FALSE}
plot(time2,temp_ensemble[1,],type='n',ylim=range(temp_ensemble),xlab="time",ylab="Temperature")
for(i in 1:NE){  
  lines(time2,temp_ensemble[i,],lwd=0.5,col="grey")
}

## sample driver rows 
drow = sample.int(nrow(precip_ensemble),Nmc,replace=TRUE)

Y.IPD <- forecastY(IC=IC[prow,"N[6,...]"],  ## sample IC
                   precip=precip_ensemble[drow,],   ## Sample drivers
                   temp=temp_ensemble[drow,],
                   beta=params[prow,"beta"],   #unsure about parameters 
                   alpha=params[prow,"alpha_site[6]"],
                   a_add=0,
                   n=Nmc)

## Plot run
plot.run()
Y.IPD.ci = apply(Y.IPD,2,quantile,c(0.025,0.5,0.975))
ecoforecastR::ciEnvelope(time2,Y.IPD.ci[1,],Y.IPD.ci[3,],col=col.alpha(N.cols[3],trans))
ecoforecastR::ciEnvelope(time2,Y.IP.ci[1,],Y.IP.ci[3,],col=col.alpha(N.cols[2],trans))
ecoforecastR::ciEnvelope(time2,Y.I.ci[1,],Y.I.ci[3,],col=col.alpha(N.cols[1],trans))
lines(time2,Y.I.ci[2,],lwd=0.5)



## Process Uncertainty 

##process error samples 
a_addmc <- 1/sqrt(params[prow,"a_add"])  ## convert from precision to standard deviation

Y.IPDE <- forecastY(IC=IC[prow,"Y[6,...]"],  ## sample IC
                    precip=precip_ensemble[drow,],   ## Sample drivers
                    temp=temp_ensemble[drow,],
                    beta=params[prow,"beta"],   #unsure about parameters 
                    alpha=params[prow,"alpha_site[6]"],
                    a_add=a_addmc,
                    n=Nmc)

## Plot run
plot.run()
Y.IPDE.ci = apply(Y.IPDE,2,quantile,c(0.025,0.5,0.975))
ecoforecastR::ciEnvelope(time2,Y.IPDE.ci[1,],Y.IPDE.ci[3,],col=col.alpha(N.cols[4],trans))
ecoforecastR::ciEnvelope(time2,Y.IPD.ci[1,],Y.IPD.ci[3,],col=col.alpha(N.cols[3],trans))
ecoforecastR::ciEnvelope(time2,Y.IP.ci[1,],Y.IP.ci[3,],col=col.alpha(N.cols[2],trans))
ecoforecastR::ciEnvelope(time2,Y.I.ci[1,],Y.I.ci[3,],col=col.alpha(N.cols[1],trans))
lines(time2,Y.I.ci[2,],lwd=0.5)








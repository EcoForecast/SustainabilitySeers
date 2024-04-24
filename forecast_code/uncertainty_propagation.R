
library(purrr)
load("~/SustainabilitySeers/data_download_code/data/ensembleParameters.Rdata")
source("~/SustainabilitySeers/Data_Download_Functions/GEFS_download.R")



#function for the forecasts.
nee_forecast <- function(ensemble) {
  params <- ensemble$params
  #grab met
  met <- ensemble$met
  temp <- met$pred_daily[which(met$variable == "air_temperature")]
  precip <- met$pred_daily[which(met$variable == "precipitation_flux")]
  humid <- met$pred_daily[which(met$variable == "relative_humidity")]
  Pres <- met$pred_daily[which(met$variable == "air_pressure")]
  LW <- met$pred_daily[which(met$variable == "surface_downwelling_shortwave_flux_in_air")]
  SW <- met$pred_daily[which(met$variable == "surface_downwelling_longwave_flux_in_air")]
  #forecast
  x_ic <- ensemble$x_ic
  mu <- rep(NA, length(temp))
  mu[1] <- x_ic
  for (t in 2:length(temp)) {
    new_nee <- mu[t-1]  + 
      params["betaX"]*mu[t-1] + 
      params["betaIntercept"] + 
      params["betaTemp"]*temp[t] + 
      params["betaPrecip"]*precip[t] + 
      params["betahumid"]*humid[t] +
      params["betaSWFlux"]*SW[t] +
      params["betaPress"]*Pres[t] +
      params["betaLWFlux"]*LW[t]
    mu[t] <- rnorm(1, new_nee, 1/sqrt(params["tau_add"]))
  }
  mu
}
#we can also setup ensemble sizes and sample from the parameter spaces.
ens <- 1000
#grab met data.
#accessing GEFS.
time_points <- seq(as.Date("2024-01-01"), as.Date("2024-01-31"), "1 day")
met_variables <- c("precipitation_flux", 
                   "air_temperature",
                   "air_pressure",
                   "relative_humidity", 
                   "surface_downwelling_shortwave_flux_in_air",
                   "surface_downwelling_longwave_flux_in_air")

site_ensemble <- vector("list", length(params))
#loop over sites
for (i in seq_along(params)) {
  #prep GEFs met.
  met <- list()
  print(paste0("Downloading GEFS weather forecasts from ", 
               time_points[1], " to ", 
               time_points[length(time_points)], 
               " for ", params[[i]]$site_id))
  #download GEFs
  pb <- utils::txtProgressBar(min = 0, max = length(time_points), style = 3)
  for (j in seq_along(time_points)) {
    met[[j]] <- GEFS_download(date = time_points[j], site_name = params[[i]]$site_id, variables = met_variables, is.daily = T)
    utils::setTxtProgressBar(pb, j)
  }
  met <- do.call(rbind, met)
  #write parameters into ensembles
  ENS <- vector("list", ens)
  for (j in seq_along(ENS)) {
    #sample met data
    ens.met <- met[which(met$parameter == sample(1:31, 1)),]
    ENS[[j]] <- list(params = params[[i]]$params[j,],
                     met = ens.met,
                     x_ic = params[[i]]$predict[j])
  }
  #forecast
  mu <- ENS %>% purrr::map(nee_forecast) %>% dplyr::bind_cols() %>% as.data.frame() %>% `colnames<-`(time_points)
  #store outputs
  site_ensemble[[i]] <- list(data = ENS, forecast = mu)
}

#time series plot.
#take BART for example
mu <- site_ensemble[[2]]$forecast
ci <- apply(mu,1,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale
plot(time_points, ci[3,], type="l", ylim=c(min(ci), max(ci)), xlab = "Date", ylab="NEE", main="NEE Forecasts")
lines(time_points, ci[1,])
lines(time_points, ci[2,], col=2)



##settings
s <- 1
Nmc = 1000  #number of Monte Carlow draws
ylim = c(min(ci), max(ci))   #not sure what to set limits as
N.cols <- c("black","red","green","blue","orange") ## set colors
trans <- 0.8       ## set transparancy
time = 1:length(temp)*2    ## total time
time1 = 1:length(temp)       ## calibration period
time2 = time1+length(temp)   ## forecast period



### Deterministic Prediction

## calculate mean of all inputs
precip.mean <- matrix(apply(precip_ensemble,2,mean),1,length(temp)) ## drivers
temp.mean <- matrix(apply(temp_ensemble,2,mean),1,length(temp))
humid.mean <- matrix(apply(humid_ensemble,2,mean),1,length(temp))
LWF.mean <- matrix(apply(LWF_ensemble,2,mean),1,length(temp))
SWF.mean <- matrix(apply(SWF_ensemble,2,mean),1,length(temp))
## parameters
#params <- do.call(rbind, ef.out$params)
params <- ensemble$params
param.mean <- apply(params,2,mean)
## initial conditions
x_ic <- ensemble$x_ic

Y.det <- nee_forecast(x_ic=mean(nee_daily, na.rm = T),  #unsure of length #do I change forecastY to forecastmu? 
                   precip=precip.mean,
                   temp = temp.mean,
                   humid = humid.mean,
                   LWF = LWF.mean,
                   SWF = SWF.mean,
                   betaX=param.mean["betaX"]*mu[t-1], 
                   betaIntercept=param.mean["betaIntercept"],
                   betaTemp=param.mean["betaTemp"]*temp[t],
                   betaPrecip=param.mean["betaPrecip"]*precip[t],
                   betahumid=param.mean["betahumid"]*humid[t],
                   betaSWFlux=param.mean["betaSWFlux"]*SWFlux[t],
                   betaPress=param.mean["betaPress"]*Press[t],
                   betaLWFlux=param.mean["betaLWFlux"]*LWFlux[t],
                   a_add=0,  #process error off
                   n=1)

## Plot run
plot.run()
lines(time2,Y.det,col="purple",lwd=3)



##Monte Carlo Propogation

##sample parameter rows from previous analysis
prow = sample.int(nrow(params),Nmc,replace=TRUE)

Y.I <-  nee_forecast(x_ic=prow(nee_daily, na.rm = T),  ## sample IC (posterier distribution, site 6 and unsure what year to pick) 
                  precip=precip.mean,
                  temp=temp.mean,
                  humid=humid.mean,
                  LWF = LWF.mean,
                  SWF = SWF.mean,
                  betaX=param.mean["betaX"]*mu[t-1], 
                  betaIntercept=param.mean["betaIntercept"],
                  betaTemp=param.mean["betaTemp"]*temp[t],
                  betaPrecip=param.mean["betaPrecip"]*precip[t],
                  betahumid= param.mean["betahumid"]*humid[t],
                  betaSWFlux=param.mean["betaSWFlux"]*SWFlux[t],
                  betaPress=param.mean["betaPress"]*Press[t],
                  betaLWFlux=param.mean["betaLWFlux"]*LWFlux[t],
                  a_add=0,  #process error off
                  n=Nmc)

  #plot run 
plot.run()
Y.I.ci = apply(Y.I,2,quantile,c(0.025,0.5,0.975))
ecoforecastR::ciEnvelope(time2,Y.I.ci[1,],Y.I.ci[3,],col=col.alpha(N.cols[1],trans))
lines(time2,Y.I.ci[2,],lwd=0.5)



##Parameter Uncertainty 

Y.IP <- nee_forecast(x_ic=prow(nee_daily, na.rm = T),  ## sample IC
                  precip=precip.mean,
                  temp=temp.mean,
                  humid=humid.mean,
                  LWF = LWF.mean,
                  SWF = SWF.mean,
                  betaX=params["betaX"]*mu[t-1], 
                  betaIntercept=params["betaIntercept"],
                  betaTemp=params["betaTemp"]*temp[t],
                  betaPrecip=params["betaPrecip"]*precip[t],
                  betahumid= params["betahumid"]*humid[t],
                  betaSWFlux=params["betaSWFlux"]*SWFlux[t],
                  betaPress=params["betaPress"]*Press[t],
                  betaLWFlux=params["betaLWFlux"]*LWFlux[t],
                  a_add=0,  #process error off
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
plot(time2,precip_ensemble[1,],type='n',ylim=range(precip_ensemble),xlab="Time",ylab="precipitation (mm)")
for(i in 1:NE){  
  lines(time2,precip_ensemble[i,],lwd=0.5,col="grey")
}


#forecast temperature 
#{r,echo=FALSE}
plot(time2,temp_ensemble[1,],type='n',ylim=range(temp_ensemble),xlab="Time",ylab="Temperature")
for(i in 1:NE){  
  lines(time2,temp_ensemble[i,],lwd=0.5,col="grey")
}

#forecast relative humidity
#{r,echo=FALSE}
plot(time2,humid_ensemble[1,],type='n',ylim=range(humid_ensemble),xlab="Time",ylab="Relative Humidity")
for(i in 1:NE){  
  lines(time2,humid_ensemble[i,],lwd=0.5,col="grey")
}

#forecast LWF
#{r,echo=FALSE}
plot(time2,LWF_ensemble[1,],type='n',ylim=range(LWf_ensemble),xlab="Time",ylab="LWFlux")
for(i in 1:NE){  
  lines(time2,LWF_ensemble[i,],lwd=0.5,col="grey")
}

#forecast SWF
#{r,echo=FALSE}
plot(time2,SWF_ensemble[1,],type='n',ylim=range(SWf_ensemble),xlab="Time",ylab="SWFlux")
for(i in 1:NE){  
  lines(time2,SWF_ensemble[i,],lwd=0.5,col="grey")
}


## sample driver rows 
drow = sample.int(nrow(precip_ensemble),Nmc,replace=TRUE)

Y.IPD <- nee_forecast(x_ic=prow(nee_daily, na.rm = T),  ## sample IC
                   precip=precip_ensemble[drow],
                   temp=temp_ensemble[drow],
                   humid=humid_ensemble[drow],
                   LWF = LWF_ensemble[drow],
                   SWF = SWF_ensemble[drow],
                   betaX=params["betaX"]*mu[t-1], 
                   betaIntercept=params["betaIntercept"],
                   betaTemp=params["betaTemp"]*temp[t],
                   betaPrecip=params["betaPrecip"]*precip[t],
                   betahumid= params["betahumid"]*humid[t],
                   betaSWFlux=params["betaSWFlux"]*SWFlux[t],
                   betaPress=params["betaPress"]*Press[t],
                   betaLWFlux=params["betaLWFlux"]*LWFlux[t],
                   a_add=0,  #process error off
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

Y.IPDE <- nee_forecast(x_ic=prow(nee_daily, na.rm = T),  ## sample IC
                    precip=precip_ensemble[drow],
                    temp=temp_ensemble[drow],
                    humid=humid_ensemble[drow],
                    LWF = LWF_ensemble[drow],
                    SWF = SWF_ensemble[drow],
                    betaX=params["betaX"]*mu[t-1], 
                    betaIntercept=params["betaIntercept"],
                    betaTemp=params["betaTemp"]*temp[t],
                    betaPrecip=params["betaPrecip"]*precip[t],
                    betahumid= params["betahumid"]*humid[t],
                    betaSWFlux=params["betaSWFlux"]*SWFlux[t],
                    betaPress=params["betaPress"]*Press[t],
                    betaLWFlux=params["betaLWFlux"]*LWFlux[t],
                    a_add=a_addmc,  #process error
                    n=Nmc)

## Plot run
plot.run()
Y.IPDE.ci = apply(Y.IPDE,2,quantile,c(0.025,0.5,0.975))
ecoforecastR::ciEnvelope(time2,Y.IPDE.ci[1,],Y.IPDE.ci[3,],col=col.alpha(N.cols[4],trans))
ecoforecastR::ciEnvelope(time2,Y.IPD.ci[1,],Y.IPD.ci[3,],col=col.alpha(N.cols[3],trans))
ecoforecastR::ciEnvelope(time2,Y.IP.ci[1,],Y.IP.ci[3,],col=col.alpha(N.cols[2],trans))
ecoforecastR::ciEnvelope(time2,Y.I.ci[1,],Y.I.ci[3,],col=col.alpha(N.cols[1],trans))
lines(time2,Y.I.ci[2,],lwd=0.5)








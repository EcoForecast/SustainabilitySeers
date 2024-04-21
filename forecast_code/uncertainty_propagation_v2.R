library(purrr)
load("~/SustainabilitySeers/data_download_code/data/ensembleParameters.Rdata")
source("~/SustainabilitySeers/Data_Download_Functions/GEFS_download.R")

#Load forecast
source("~/SustainabilitySeers/forecast_code/03_forecast_v2.R")

  

##settings
s <- 1
Nmc = 1000  #number of Monte Carlow draws
ylim = c(50,1000)   #not sure what to set limits as
N.cols <- c("black","red","green","blue","orange") ## set colors
trans <- 0.8       ## set transparancy
time = 1:length(temp)*2    ## total time
time1 = 1:length(temp)       ## calibration period
time2 = time1+length(temp)   ## forecast period

plot.run <- function(){
  sel = seq(s,ncol(ci),by=NS)
  plot(time,time,type='n',ylim=ylim,ylab="N")
  ecoforecastR::ciEnvelope(time1,ci[1,sel],ci[3,sel],col=col.alpha("lightBlue",0.6))
  lines(time1,ci[2,sel],col="blue")
  points(time1,No[s,])
}

#ci <- apply(as.matrix(out$predict),2,quantile,c(0.025,0.5,0.975))



### Deterministic Prediction

## calculate mean of all inputs
precip.mean <- matrix(apply(precip,2,mean),1,length(temp)) ## drivers
temp.mean <- matrix(apply(temp,2,mean),1,length(temp))
humid.mean <- matrix(apply(humid,2,mean),1,length(temp))
LW.mean <- matrix(apply(LW,2,mean),1,length(temp))
SW.mean <- matrix(apply(SW,2,mean),1,length(temp))
Press.mean <- matrix(apply(Press,2,mean),1,length(temp))
## parameters
params <- ensemble$params
param.mean <- apply(params,2,mean)
## initial conditions
x_ic <- ensemble$x_ic

Y.det <- nee_forecast(x_ic=mean(x_ic[,"mu[1,30"]),  #unsure of length 
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
                      tau_add=0,  #process error off
                      n=1)

## Plot run
plot.run()
lines(time2,Y.det,col="purple",lwd=3)



##Monte Carlo Propogation

##sample parameter rows from previous analysis
prow = sample.int(nrow(params),Nmc,replace=TRUE)

Y.I <-  nee_forecast(x_ic=x_ic[prow,"mu[1,30]"],  ## sample IC (posterier distribution, site 1 and unsure what year to pick) 
                     precip=precip.mean,
                     temp=temp.mean,
                     humid=humid.mean,
                     LW = LW.mean,
                     SW = SW.mean,
                     Press = Press.mean,
                     betaX=param.mean["betaX"]*mu[t-1], 
                     betaIntercept=param.mean["betaIntercept"],
                     betaTemp=param.mean["betaTemp"]*temp[t],
                     betaPrecip=param.mean["betaPrecip"]*precip[t],
                     betahumid= param.mean["betahumid"]*humid[t],
                     betaSWFlux=param.mean["betaSWFlux"]*SWFlux[t],
                     betaPress=param.mean["betaPress"]*Press[t],
                     betaLWFlux=param.mean["betaLWFlux"]*LWFlux[t],
                     tau_add=0,  #process error off
                     n=Nmc)

#plot run 
plot.run()
Y.I.ci = apply(Y.I,2,quantile,c(0.025,0.5,0.975))
ecoforecastR::ciEnvelope(time2,Y.I.ci[1,],Y.I.ci[3,],col=col.alpha(N.cols[1],trans))
lines(time2,Y.I.ci[2,],lwd=0.5)




##Parameter Uncertainty 

Y.IP <- nee_forecast(x_ic=x_ic[prow, "mu[1,30]"],  ## sample IC
                     precip=precip.mean,
                     temp=temp.mean,
                     humid=humid.mean,
                     LW = LW.mean,
                     SW = SW.mean,
                     Press = Press.mean,
                     betaX=params[prow, "betaX"], 
                     betaIntercept=params[prow,"betaIntercept"],
                     betaTemp=params[prow, "betaTemp"],
                     betaPrecip=params[prow,"betaPrecip"],
                     betahumid= params[prow,"betahumid"],
                     betaSWFlux=params[prow,"betaSWFlux"],
                     betaPress=params[prow,"betaPress"],
                     betaLWFlux=params[prow,"betaLWFlux"],
                     tau_add=0,  #process error off
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
plot(time2,precip[1,],type='n',ylim=range(precip),xlab="Time",ylab="precipitation (mm)")
for(i in 1:NE){  
  lines(time2,precip[i,],lwd=0.5,col="grey")
}


#forecast temperature 
#{r,echo=FALSE}
plot(time2,temp[1,],type='n',ylim=range(temp),xlab="Time",ylab="Temperature")
for(i in 1:NE){  
  lines(time2,temp[i,],lwd=0.5,col="grey")
}

#forecast relative humidity
#{r,echo=FALSE}
plot(time2,humid[1,],type='n',ylim=range(humid),xlab="Time",ylab="Relative Humidity")
for(i in 1:NE){  
  lines(time2,humid[i,],lwd=0.5,col="grey")
}

#forecast LW
#{r,echo=FALSE}
plot(time2,LW[1,],type='n',ylim=range(LW),xlab="Time",ylab="LWFlux")
for(i in 1:NE){  
  lines(time2,LW[i,],lwd=0.5,col="grey")
}

#forecast SW
#{r,echo=FALSE}
plot(time2,SW[1,],type='n',ylim=range(SW),xlab="Time",ylab="SWFlux")
for(i in 1:NE){  
  lines(time2,SW[i,],lwd=0.5,col="grey")
}

#forecast air pressure
#{r,echo=FALSE}
plot(time2,Press[1,],type='n',ylim=range(Press),xlab="Time",ylab="SWFlux")
for(i in 1:NE){  
  lines(time2,Press[i,],lwd=0.5,col="grey")
}

## sample driver rows 
drow = sample.int(nrow(precip),Nmc,replace=TRUE)

Y.IPD <- nee_forecast(x_ic=x_ic[prow, "mu[1,30]"],  ## sample IC
                      precip=precip[drow],
                      temp=temp[drow],
                      humid=humid[drow],
                      LW = LW[drow],
                      SW = SW[drow],
                      Press = Press[drow],
                      betaX=params[prow, "betaX"], 
                      betaIntercept=params[prow,"betaIntercept"],
                      betaTemp=params[prow, "betaTemp"],
                      betaPrecip=params[prow,"betaPrecip"],
                      betahumid= params[prow,"betahumid"],
                      betaSWFlux=params[prow,"betaSWFlux"],
                      betaPress=params[prow,"betaPress"],
                      betaLWFlux=params[prow,"betaLWFlux"],
                      tau_add=0,  #process error off
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
tau_addmc <- 1/sqrt(params[prow,"tau_add"])  ## convert from precision to standard deviation

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
                       tau_add=tau_addmc,  #process error
                       n=Nmc)

## Plot run
plot.run()
Y.IPDE.ci = apply(Y.IPDE,2,quantile,c(0.025,0.5,0.975))
ecoforecastR::ciEnvelope(time2,Y.IPDE.ci[1,],Y.IPDE.ci[3,],col=col.alpha(N.cols[4],trans))
ecoforecastR::ciEnvelope(time2,Y.IPD.ci[1,],Y.IPD.ci[3,],col=col.alpha(N.cols[3],trans))
ecoforecastR::ciEnvelope(time2,Y.IP.ci[1,],Y.IP.ci[3,],col=col.alpha(N.cols[2],trans))
ecoforecastR::ciEnvelope(time2,Y.I.ci[1,],Y.I.ci[3,],col=col.alpha(N.cols[1],trans))
lines(time2,Y.I.ci[2,],lwd=0.5)
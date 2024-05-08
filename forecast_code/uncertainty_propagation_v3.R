library(purrr)
load("~/SustainabilitySeers/data_download_code/data/ensembleParameters.Rdata")
source("~/SustainabilitySeers/Data_Download_Functions/GEFS_download.R")

#Load forecast
source("~/SustainabilitySeers/forecast_code/forecast_function.R")
source("~/SustainabilitySeers/forecast_code/03_forecast_v2.R")
source("~/SustainabilitySeers/calibration_code/02_DLM_Calibration.R")


1##settings
s <- 1
Nmc = 1000  #number of Monte Carlow draws
ylim = c(50,1000)   #not sure what to set limits as
N.cols <- c("black","red","green","blue","orange") ## set colors
trans <- 0.8       ## set transparancy
time = 1:time.bart*2    ## total time
time1 = 1:time.bart       ## calibration period
time2 = time1+time.bart   ## forecast period

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
ap.mean  <- rep(mean(ap.bart), length(bart))
temp.mean <- rep(mean(temp.bart), length(bart))
precip.mean <- rep(mean(precip.bart), length(bart))
humid.mean <- rep(mean(humid.bart), length(bart))
LWF.mean <- rep(mean(LWF.bart), length(bart))
SWF.mean <- rep(mean(SWF.bart), length(bart))

#precip.mean <- matrix(apply(precip,2,mean),1,length(temp)) ## drivers
#temp.mean <- matrix(apply(temp,2,mean),1,length(temp))
#humid.mean <- matrix(apply(humid,2,mean),1,length(temp))
#LW.mean <- matrix(apply(LW,2,mean),1,length(temp))
#SW.mean <- matrix(apply(SW,2,mean),1,length(temp))
#Press.mean <- matrix(apply(Press,2,mean),1,length(temp))
## parameters
#params <- ensemble$params
#param.mean <- apply(params,2,mean)
## initial conditions
x_ic = mean(nee_daily.bart, na.rm = T)
#x_ic <- ensemble$x_ic

Y.det <- list(x_ic = mean(nee_daily.bart, na.rm = T),
              tau_ic = 0,  # uncertainty of initial condition off
              a_obs = 1,
              r_obs = 1,
              a_add = 1,
              r_add = 0,  #process error off 
              n = length(nee_daily.bart),
              y = nee_daily.bart,
              Press = ap.mean,
              Precip = precip.mean,
              Temp = temp.mean,
              humid = humid.mean,
              LWFlux = LWF.mean,
              SWFlux = SWF.mean)

## Plot run
plot.run()
lines(time2,Y.det,col="purple",lwd=3)



##Monte Carlo Propogation

##sample parameter rows from previous analysis
prow = sample.int(nrow(params),Nmc,replace=TRUE)

Y.I <-  list(x_ic = mean(nee_daily.bart, na.rm = T),
             tau_ic = 1/sd(nee_daily.bart, na.rm = T),   #IC uncertainty on 
             a_obs = 1,
             r_obs = 1,
             a_add = 1,
             r_add = 0,  #process error off 
             n = length(nee_daily.bart),
             y = nee_daily.bart,
             Press = ap.mean,
             Precip = precip.mean,
             Temp = temp.mean,
             humid = humid.mean,
             LWFlux = LWF.mean,
             SWFlux = SWF.mean)

#plot run 
plot.run()
Y.I.ci = apply(Y.I,2,quantile,c(0.025,0.5,0.975))
ecoforecastR::ciEnvelope(time2,Y.I.ci[1,],Y.I.ci[3,],col=col.alpha(N.cols[1],trans))
lines(time2,Y.I.ci[2,],lwd=0.5)




##Parameter Uncertainty 

Y.IP <- list(x_ic = mean(nee_daily.bart, na.rm = T),
             tau_ic = 1/sd(nee_daily.bart, na.rm = T),
             a_obs = 1,
             r_obs = 1,
             a_add = 1,
             r_add = 0,
             n = length(nee_daily.bart),
             y = nee_daily.bart,
             Press = ap.mean,
             Precip = precip.mean,
             Temp = temp.mean,
             humid = humid.mean,
             LWFlux = LWF.mean,
             SWFlux = SWF.mean)

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
plot(time2,precip.bart[1,],type='n',ylim=range(precip.bart),xlab="Time",ylab="precipitation (mm)")
for(i in 1:NE){  
  lines(time2,precip.bart[i,],lwd=0.5,col="grey")
}


#forecast temperature 
#{r,echo=FALSE}
plot(time2,temp.bart[1,],type='n',ylim=range(temp.bart),xlab="Time",ylab="Temperature")
for(i in 1:NE){  
  lines(time2,temp.bart[i,],lwd=0.5,col="grey")
}

#forecast relative humidity
#{r,echo=FALSE}
plot(time2,humid.bart[1,],type='n',ylim=range(humid.bart),xlab="Time",ylab="Relative Humidity")
for(i in 1:NE){  
  lines(time2,humid.bart[i,],lwd=0.5,col="grey")
}

#forecast LW
#{r,echo=FALSE}
plot(time2,LWF.bart[1,],type='n',ylim=range(LWF.bart),xlab="Time",ylab="LWFlux")
for(i in 1:NE){  
  lines(time2,LWF.bart[i,],lwd=0.5,col="grey")
}

#forecast SW
#{r,echo=FALSE}
plot(time2,SWF.bart[1,],type='n',ylim=range(SWF.bart),xlab="Time",ylab="SWFlux")
for(i in 1:NE){  
  lines(time2,SWF.bart[i,],lwd=0.5,col="grey")
}

#forecast air pressure
#{r,echo=FALSE}
plot(time2,ap.bart[1,],type='n',ylim=range(ap.bart),xlab="Time",ylab="SWFlux")
for(i in 1:NE){  
  lines(time2,ap.bart[i,],lwd=0.5,col="grey")
}

## sample driver rows 
drow = sample.int(nrow(precip.bart),Nmc,replace=TRUE)

Y.IPD <- list(x_ic = mean(nee_daily.bart, na.rm = T),
              tau_ic = 1/sd(nee_daily.bart, na.rm = T),
              a_obs = 1,
              r_obs = 1,
              a_add = 1,
              r_add = 0,
              n = length(nee_daily.bart),
              y = nee_daily.bart,
              Press = ap.bart,
              Precip = precip.bart,
              Temp = temp.bart,
              humid = humid.bart,
              LWFlux = LWF.bart,
              SWFlux = SWF.bart)

## Plot run
plot.run()
Y.IPD.ci = apply(Y.IPD,2,quantile,c(0.025,0.5,0.975))
ecoforecastR::ciEnvelope(time2,Y.IPD.ci[1,],Y.IPD.ci[3,],col=col.alpha(N.cols[3],trans))
ecoforecastR::ciEnvelope(time2,Y.IP.ci[1,],Y.IP.ci[3,],col=col.alpha(N.cols[2],trans))
ecoforecastR::ciEnvelope(time2,Y.I.ci[1,],Y.I.ci[3,],col=col.alpha(N.cols[1],trans))
lines(time2,Y.I.ci[2,],lwd=0.5)



## Process Uncertainty 

##process error samples 
#tau_addmc <- 1/sqrt(params[prow,"tau_add"])  ## convert from precision to standard deviation

Y.IPDE <- list(x_ic = mean(nee_daily.bart, na.rm = T),
               tau_ic = 1/sd(nee_daily.bart, na.rm = T),
               a_obs = 1,
               r_obs = 1,
               a_add = 1,
               r_add = 1,   #process error on 
               n = length(nee_daily.bart),
               y = nee_daily.bart,
               Press = ap.bart,
               Precip = precip.bart,
               Temp = temp.bart,
               humid = humid.bart,
               LWFlux = LWF.bart,
               SWFlux = SWF.bart)

## Plot run
plot.run()
Y.IPDE.ci = apply(Y.IPDE,2,quantile,c(0.025,0.5,0.975))
ecoforecastR::ciEnvelope(time2,Y.IPDE.ci[1,],Y.IPDE.ci[3,],col=col.alpha(N.cols[4],trans))
ecoforecastR::ciEnvelope(time2,Y.IPD.ci[1,],Y.IPD.ci[3,],col=col.alpha(N.cols[3],trans))
ecoforecastR::ciEnvelope(time2,Y.IP.ci[1,],Y.IP.ci[3,],col=col.alpha(N.cols[2],trans))
ecoforecastR::ciEnvelope(time2,Y.I.ci[1,],Y.I.ci[3,],col=col.alpha(N.cols[1],trans))
lines(time2,Y.I.ci[2,],lwd=0.5)

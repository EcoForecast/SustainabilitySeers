# PREAMBLE
# Ecological Forecasting Milestone 5:
# Historical Fit
# Date Updated: ''

# Load Libraries
library(neon4cast)
library(dplyr)
library(lubridate)
library(ggplot2)
library(arrow)
library(glue)
library(readr)
library(ecoforecastR)
library(rjags)


#Get the data ----
source("~/SustainabilitySeers/Data_Download_Functions/01_datatargetdownload.R")

met_future <- met_future %>%
  rename(prediction=observation)

met_all <- rbind(met_future, df_past)
nee <- target1 %>% filter(variable=="nee") 
le <- target1 %>% filter(variable=="le")

# Look at quick time series ----
plot(x=nee$datetime, y=nee$observation,type='l',ylab="NEE",lwd=2)
plot(x=le$datetime, y=le$observation,type='l',ylab="LE",lwd=2, col=2)


# Create a random walk model ----
## Define the model ----
RandomWalk = "
model{
  
  #### Data Model
  for(t in 1:n){
    y[t] ~ dnorm(x[t],tau_obs)
  }
  
  #### Process Model
  for(t in 2:n){
    x[t]~dnorm(x[t-1],tau_add)
  }
  
  #### Priors
  x[1] ~ dnorm(x_ic,tau_ic)
  tau_obs ~ dgamma(a_obs,r_obs)
  tau_add ~ dgamma(a_add,r_add)
}
"

# Define the data (NEE) ----
time <- nee$datetime
y <- nee$observation
hist(y) # this is mostly normal although high variance, don't think we need a log transform

data <- list(y=y,n=length(y),      ## data
             x_ic=log(1000),tau_ic=100, ## initial condition prior
             a_obs=1,r_obs=1,           ## obs error prior
             a_add=1,r_add=1            ## process error prior
)

nchain = 3
init <- list()
for(i in 1:nchain){
  y.samp = sample(y,length(y),replace=TRUE)
  init[[i]] <- list(tau_add=1/var(diff(y.samp)),  ## initial guess on process precision
                    tau_obs=5/var(y.samp))        ## initial guess on obs precision
}


# Now run the model ----
gc()
nee.model   <- jags.model (file = textConnection(RandomWalk),
                         data = data,
                         inits = init,
                         n.chains = 3)
n.out   <- coda.samples (model = nee.model,
                            variable.names = c("x","tau_add","tau_obs"),
                            n.iter = 10000)

# Run diagnostics ----
plot(n.out) # check for convergence
effectiveSize(n.out)

## check variable correlations 
pairs(n.out)
cor(n.out)

summary(n.out)

# Now plot a time series ----
time.rng = c(1,length(time))       ## adjust to zoom in and out
out <- as.matrix(n.out)         ## convert from coda to matrix  
x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(out[,x.cols],2,quantile,c(0.025,0.5,0.975)) 

plot(time,ci[2,],type='n',ylim=range(y,na.rm=TRUE),ylab="NEE",xlim=time[time.rng])

## adjust x-axis label to be monthly if zoomed
if(diff(time.rng) < 100){ 
  axis.Date(1, at=seq(time[time.rng[1]],time[time.rng[2]],by='month'), format = "%Y-%m")
}
ecoforecastR::ciEnvelope(time,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(time,y,pch="+",cex=0.5)










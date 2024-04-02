load("~/SustainabilitySeers/Data/calibration.Rdata")
#set initial value of NEE for the forecasts.
x_ic <- 0.5177017
#combine parameters from ef.out object.
params <- do.call(rbind, ef.out$params)
#we can also setup ensemble sizes and sample from the parameter spaces.
if (FALSE) {
  ens <- 100
  sample.inds <- sample(seq_along(params[,1]), ens)
  params <- params[sample.inds,]
}
#grab met data.
met <- as.data.frame(ef.out$data$Xf)
#remove na.
na.inds <- which(is.na(met[,2]))
met <- met[-na.inds,]
#function for the forecasts.
nee_forecast <- function(params, met, x_ic) {
  mu <- matrix(NA, dim(params)[1], dim(met)[1])
  mu[,1] <- rep(x_ic, dim(params)[1])
  for (t in 2:dim(met)[1]) {
    new_nee <- mu[,t-1]  + 
      params[,"betaX"]*mu[,t-1] + 
      params[,"betaIntercept"] + 
      params[,"betaTemp"]*met$Temp[t] + 
      params[,"betaPrecip"]*met$Precip[t] + 
      params[,"betahumid"]*met$humid[t]
    mu[,t] <- new_nee
  }
  mu
}
#run forecasts.
mu <- nee_forecast(params, met, x_ic)
#time series plot.
ci <- apply(mu,2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale
plot(ci[3,], type="l", ylim=c(min(ci), max(ci)), xlab = "DOY", ylab="NEE", main="NEE Forecasts")
lines(ci[1,])
lines(ci[2,], col=2)

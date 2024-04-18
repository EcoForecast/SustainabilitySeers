# Define forecast function


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
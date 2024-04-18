#load ensemble parameters
load("~/SustainabilitySeers/data_download_code/data/ensembleParameters.Rdata")
source("~/SustainabilitySeers/Data_Download_Functions/GEFS_download.R")

# Load covariate data ----
source("~/SustainabilitySeers/data_download_code/01_datatargetdownload.R") # NEE
source("~/SustainabilitySeers/data_download_code/01A_NOAA_datadownload.R") # weather

#Load forecast function
source("~/SustainabilitySeers/forecast_code/forecast_function.R")

# Load forcast outputs
source("~/SustainabilitySeers/forecast_code/03_forecast_v2.R")

##' Kalman Filter: Analysis step
##' @param  mu.f = Forecast mean (vector)
##' @param  P.f  = Forecast covariance (matrix)
##' @param  Y    = observations, with missing values as NAs) (vector)
##' @param  R    = observation error covariance (matrix)
##' @param  H    = observation matrix (maps observations to states)
KalmanAnalysis <- function(mu.f,P.f,Y,R,H,I){
  obs = !is.na(Y) ## which Y's were observed?
  if(any(obs)){
    H <- H[obs,]                                              ## observation matrix
    if (dim(t(H))[1] == 1) {
      K <- P.f %*% t(t(H)) %*% solve(H%*%P.f%*%t(t(H)) + R[obs,obs])  ## Kalman gain
    } else{
      K <- P.f %*% t(H) %*% solve(H%*%P.f%*%t(H) + R[obs,obs])  ## Kalman gain
    }
    mu.a <- mu.f + K%*%(Y[obs] - H %*% mu.f)                  ## update mean
    P.a <- (I - K %*% H)%*%P.f                                ## update covariance
    ## Note: Here's an alternative form that doesn't use the Kalman gain
    ## it is less efficient due to the larger number of matrix inversions (i.e. solve)
    ## P.a <- solve(t(H)%*%solve(R[obs,obs])%*%(H) + solve(P.f))                             
    ## mu.a <- P.a %*% (t(H)%*%solve(R[obs,obs])%*%Y[obs] + solve(P.f)%*%mu.f)
  } else {
    ##if there's no data, the posterior is the prior
    mu.a = mu.f
    P.a = P.f
  }
  return(list(mu.a=mu.a,P.a=P.a))
}

# Start data assimilation
# define period
nt <- length(time_points)
nee <- target1 %>% filter(variable=="nee")
nsite <- length(site)
site_ids <- site

# filter nee observations
nee.filter <- nee %>% dplyr::filter(datetime %in% time_points) 

# convert nee to a more generalized data frame
nee.mean <- nee.sd <- data.frame(matrix(NA, nrow = nsite, ncol = nt)) %>%
  `rownames<-`(site) %>%
  `colnames<-`(time_points)
## assume constant standard error for test run. ##
nee.sd[,] <- 0.2
for (id in seq_along(site)) {
  ind <- which(nee.filter$datetime %in% time_points & nee.filter$site_id == site[id])
  nee.mean[id, as.character(nee.filter$datetime[ind])] <- nee.filter$observation[ind]
}

# data assimilation
ens <- 1000
mu.f  = matrix(NA,nsite,nt+1)  ## forecast mean for time t
mu.a  = matrix(NA,nsite,nt)  ## analysis mean for time t
P.f  = array(NA,c(nsite,nsite,nt+1))  ## forecast variance for time t
P.a  = array(NA,c(nsite,nsite,nt))  ## analysis variance for time t
x_ic <- site_ensemble %>% 
  purrr::map('data') %>%
  unlist(recursive = F) %>%
  purrr::map('x_ic') %>%
  unlist %>%
  matrix(nrow = ens, ncol = nsite)
mu.f[,1] <- colMeans(x_ic)
P.f[,,1] <- cov(x_ic)
for (t in seq_along(time_points)) {
  # run forecast
  x_ic <- site_ensemble %>% 
    purrr::map('data') %>%
    unlist(recursive = F) %>%
    purrr::map2(t, nee_forecast) %>%
    unlist %>%
    matrix(nrow = ens, ncol = nsite) 
  
  # update mu.f and P.f
  mu.f[,t+1] <- colMeans(x_ic)
  P.f[,,t+1] <- cov(x_ic)
  Y <- nee.mean[,t]
  R <- diag(nee.sd[,t], nsite)
  H <- I <- diag(1, nsite)
  
  # run analysis
  analysis <- KalmanAnalysis(mu.f[,t+1], P.f[,,t+1], Y, R, H, I)
  mu.a[,t] <- analysis$mu.a
  P.a[,,t] <- analysis$P.a
  
  # update site_ensemble
  update <- mvtnorm::rmvnorm(ens, analysis$mu.a, analysis$P.a, method = "svd")
  
  # update initial conditions.
  for (i in seq_along(site_ids)) {
    for (j in 1:ens) {
      site_ensemble[[i]]$data[[j]]$x_ic <- update[j, i]
    }
  }
}

# plot
for(i in 1:nsite){
  ci = rbind(mu.a[i,]-1.96*sqrt(P.a[i,i,]),mu.a[i,]+1.96*sqrt(P.a[i,i,]))
  plot(time_points,mu.a[i,],ylim=range(ci,na.rm=TRUE),type='n',main=site_ids[i],xlab="Date",ylab="NEE")
  ecoforecastR::ciEnvelope(time_points,ci[1,],ci[2,],col="lightBlue")
  lines(time_points,mu.a[i,],col=4)
  points(time_points,nee.mean[i,],col="red")
  legend("topleft", lty=c(NA,1), pch=c("o", NA), col=c("red", "lightBlue"), legend = c("Observation", "Forecast"))
}



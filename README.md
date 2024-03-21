# SustainabilitySeers

## Contact information

Dongchen Zhang
Phone Number: 617-320-9561
Email: zhangdc@bu.edu

Breanna van Loenen
bvanloen@bu.edu
3474176834

Tessa Keeney
Email: tkeeney@bu.edu 
Cell: 360-854-8174

Katherine Losada
Email: klosada@bu.edu
Phone number: (603)707-6791


#Milestone 5: SustainabilitySeers
#Dongchen, Tessa, Breanna


#Data Model
In the ef.out model, our data model relates observations (y and covariates) to the latent states of 1) NEE observations with prescribed priors on the observation error and 2) covariates (temperature, precipitation, and relative humidity) and intercept with observation errors estimated by fixed effects.

#Process Model 
The process model relates the latent state of NEE to the prior of NEE, intercept, and covariates with prescribed process error. We have chosen to model Gaussian process error, although future iterations of our model will attempt to apply Laplace distribution to account for the non-normal behavior of the NEE variable. 

#Priors 
The priors consist of the observation error and the process error, as well as the initial conditions of the model. We chose gamma distributions for our priors since they are conjugates of the normally distributed mu, but will likely alter this in the model in future iterations to reflect the Laplace distribution of NEE. Fixed effects within the priors include the means and precisions associated with the covariates (temperature, precipitation, and relative humidity). 


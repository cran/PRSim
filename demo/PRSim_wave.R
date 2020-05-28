# short demo

### for empirical distribution
out <- prsim.wave(data=runoff_multi_sites, number_sim=1, marginal="empirical")
### for kappa distribution
out <- prsim.wave(data=runoff_multi_sites, number_sim=1, marginal="kappa", GoFtest = "KS")
# out <- prsim.wave(data=runoff_multi_sites, number_sim=5, marginal="kappa",
#                   GoFtest = NULL,pars=NULL, p_val=NULL)
# simulations_multi_sites<-list(out[[1]]$simulation,out[[2]]$simulation,out[[3]]$simulation,out[[4]]$simulation)
# setwd("C:/Users/mbrunner/Documents/PRSim-devel/data")
# save(simulations_multi_sites,file='simulations_multi_sites.rda')

### for GEV distribution
require("evd")
require("ismev")
rGEV <- function(n, theta)  rgev(n, theta[1], theta[2], theta[3])
pGEV <- function(x, theta)  pgev(x, theta[1], theta[2], theta[3])
GEV_fit <- function( xdat, ...)   gev.fit( xdat, show=FALSE, ...)$mle


### GEV
out <- prsim.wave(data=runoff_multi_sites, number_sim=1, marginal="GEV", GoFtest = "KS", n_par=3)


### load stochastically simulated time series
### station 2
sim <- out[[2]]$simulation

par(mai=c(.9,.9,.1,.1))
### observed time series
plot(sim$timestamp[1:1000], sim$Qobs[1:1000], type="l", 
     xlab="Time [d]", ylab=expression(paste("Discharge [m"^3,"/s]")))
### add simulations
matlines(sim$timestamp[1:1000], sim[1:1000, grep("r", names(sim))],
         lty=1, col="gray")



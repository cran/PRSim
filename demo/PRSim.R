# short demo
data(runoff)

out <- prsim(runoff, 1, 10)
sim <- out$simulation
par(mai=c(.9,.9,.1,.1))
plot(sim$timestamp[1:1000], sim$Qobs[1:1000], type="l", 
     xlab="Time [d]", ylab=expression(paste("Discharge [m"^3,"/s]")))
matlines(sim$timestamp[1:1000], sim[1:1000, grep("r", names(sim))],
         lty=1, col="gray")



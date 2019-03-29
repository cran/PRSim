# some simple testing commands...
#   testthat would be an overshoot...

require(PRSim)

# raw demos
demo( "PRSim", ask=FALSE)
demo( "PRSim-validate", ask=FALSE)



# testing input
data(runoff)
unique(runoff$YYYY)


try( prsim( runoff[1:130,] ))   #  At least one year of data required.
try( prsim( runoff[1:730,] ))   #  At least one year of data required.
try( prsim( runoff[1:1445,] ))  #  No missing values allowed. Some days are missing.
try( prsim( runoff[runoff$YYYY<1976,] )) # At least one year of data required.



suppressWarnings( out <- prsim( runoff[runoff$YYYY<1977,] ) )

runof <- runoff[runoff$YYYY<1980,]

set.seed(1)
out1 <- prsim( runof, kappapar=FALSE, suppWarn=TRUE)

runo <- runof
names( runo) <- tolower( names(runof))
try( prsim( runo, kappapar=FALSE, suppWarn=TRUE)) #  Wrong column for observations selected.

runo <- runof[,4:1]
set.seed(1)
out3 <- prsim( runo, kappapar=FALSE, suppWarn=TRUE) #  ok
identical(out1,out3) 

runo <- runof[,4:1]
set.seed(1)
out4 <- prsim( runo, station_id=1, kappapar=FALSE, suppWarn=TRUE) #  ok
identical(out1,out4) 

tmp <- paste(runof$YYYY, runof$MM, runof$DD,sep=" ")
runo <- data.frame(time=as.POSIXct(strptime(tmp, format="%Y %m %d", tz="GMT")), Qobs=runof$Qobs)
set.seed(1)
out5 <- prsim( runo, kappapar=FALSE, suppWarn=TRUE) #  ok
identical(out1,out5) 

# 
\name{PRSim-package}
\alias{PRSim-package}
\alias{PRSim}
\docType{package}
\title{
\packageTitle{PRSim}
}
\description{
\packageDescription{PRSim}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{PRSim}
\packageIndices{PRSim}

Simulation in the frequency domain is based on the randomization of the phases of the Fourier transform. We here combine phase randomization simulation with the flexible, four-parameter Kappa distribution, which allows for the extrapolation to yet unobserved low and high flows. The
simulation approach consists of eight steps: 1) fitting of theoretical Kappa distribution, 2) normalization and deseasonalization, 3) Fourier transformation, 4) Fourier phases computation, 5) random phase generation, 6) inverse Fourier transformation, 7) back transformation, and 8) simulation.
}
\author{
\packageAuthor{PRSim}

Maintainer: \packageMaintainer{PRSim}
}
\references{
Brunner, Bardossy, Furrer (2019) Technical note: Stochastic simulation of streamflow time series using phase randomization. Submitted. 
}
\keyword{ package }
\examples{
\dontrun{
demo("PRSim")
demo("PRSim-validate")
}
}

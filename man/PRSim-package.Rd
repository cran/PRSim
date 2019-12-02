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
Contains two functions for the stochastic simulation of continuous discharge time series: prsim and prsim.wave both using phase randomization. prsim is based on the Fourier transform while prsim.wave uses the wavelet transform.

prsim: Simulation in the frequency domain is based on the randomization of the phases of the Fourier transform. We here combine phase randomization simulation with the flexible, four-parameter kappa distribution, which allows for the extrapolation to yet unobserved low and high flows. Alternative distributions or the empirical distribution can be used instead. The
simulation approach consists of eight steps: (1) fitting of theoretical Kappa distribution, (2) normalization and deseasonalization, (3) Fourier transformation, (4) Fourier phases computation, (5) random phase generation, (6) inverse Fourier transformation, (7) back transformation, and (8) simulation.

prsim.wave: Simulation in the frequency domain based on the randomization of the phases of the continuous wavelet transform. We combine phase randomization with the flexible, four-parameter kappa distribution. Alternative theoretical distributions or the empirical distribution can be used instead. The simulation procedure consists of five steps: (1) Derivation of random phases from a white noise time series, (2) Fitting of kappa distribution, (3) Wavelet transform, (4) Inverse wavelet transform, and (5) Transformation to the kappa distribution (or the distribution of choice).
}
\author{
\packageAuthor{PRSim}

Maintainer: \packageMaintainer{PRSim}
}
\references{
Brunner, M. I., A. Bárdossy, and R. Furrer (2019). Technical note: Stochastic simulation of streamflow time series using phase randomization. Hydrology and Earth System Sciences, 23, 3175-3187, doi:10.5194/hess-23-3175-2019
}
\keyword{ package }
\examples{
\dontrun{
demo("PRSim")
demo("PRSim-validate")
demo("PRSim_wave")
demo("PRSim_wave-validate")
}
}

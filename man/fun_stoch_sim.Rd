\name{pRsim}
\alias{PRsim}
\alias{prsim}
\title{Simulate for one station}
\description{
Applies the algorithm to a single station
}
\usage{
prsim(data, station_id="Qobs", number_sim=1, win_h_length=15,
      verbose=TRUE, kappapar=TRUE, suppWarn=FALSE, KStest=FALSE)
}
\arguments{
  \item{data}{data frame containing the time indications and runoff of at least one station. See \sQuote{Details}.}
  \item{station_id}{identifies the station in case several runoffs are present in \code{data}. See \sQuote{Details}.}
  \item{win_h_length}{(half-)length of moving window size.}
  \item{number_sim}{number of simulations to be carried out.}
  \item{verbose}{logical. Should progress be reported?}
  \item{kappapar}{logical. Should the estimated parameters of the kappa distribution be returned?}
  \item{suppWarn}{logical. See \sQuote{Details}.}
  \item{KStest}{logical. Should a \code{ks.test} for daily data be performed and the p-values be passed back?}
}
\details{
Time can be given with three columns named \code{"YYYY"}, \code{"MM"}, \code{"DD"}, or as in POSIXct format \code{YYYY-MM-DD}.
All leap days (Feb 29th) will be omitted from the analysis, but no missing observations are allowed.

Stations are identified by column name (default \code{"Qobs"}), or by column index.


The function \code{homtest::par.kappa} might issue quite a few warnings of type \code{In fn(par, ...) : value out of range in 'gammafn'}. The argument \code{suppWarn} allows to silence warnings for the specific function call via \code{suppressWarnings()}. Of course, a subsequent check via \code{warnings()} is recommended.


}
\value{A list with elements 
 \item{simulation }{A data frame with time information, observations, deseaonalized observations and 
       \code{number_sim} columns containing the simulated runoff.}
  \item{kappa_pars}{A matrix containing the estimated parameters of the kappa distribution (if \code{kappapar}).}
 \item{p_val}{A vector containing the p-values of \code{ks.test} applied to the daily detrended data (if \code{KStest}).}
}
\references{
Brunner, Bardossy, Furrer (2019) Technical note: Stochastic simulation of streamflow time series using phase randomization. Submitted.
}
\author{
Manuela Brunner
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data( runoff)
out <- prsim( runoff[ runoff$YYYY<1980, ], "Qobs", 1, suppWarn=TRUE)
# warnings() # as a follow-up to `suppWarn=TRUE`

}
\keyword{ts}

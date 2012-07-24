\name{coef.bayesCox}
\alias{coef.bayesCox}
\title{Extract Coefficients from Bayesian Cox Model}
\description{
  Extract coefficient values from \code{bayesCox} fitting
  results, and summarize the posterior mean, posterior 2.5\% and 97.5\%
  quantiles into a data frame.
}
\usage{
\method{coef}{bayesCox}(object, \dots)
}
\arguments{
  \item{object}{
    an object returned by function \code{bayesCox}.
  }
  \item{\dots}{
    other arguments.
  }
}
\value{
  A data.frame with 6 columns \code{("Low", "Mid", "High", "Time",
    "Cov", "Model")}, where \code{"Low"} and \code{"High"} are the
  posterior 2.5\% and 97.5\% quantiles; \code{"Mid"} is the posterior
  mean; \code{"Cov"} and \code{"Model"} contain character values of the
  covariates and model types.
}
\seealso{
  \code{\link{bayesCox}}, and \code{\link{plotCoef}}.
}
\examples{
# See the examples in bayesCox
}
\keyword{extract}
\keyword{bayesCox}
\keyword{coefficient}

\name{coef.tvTran}
\alias{coef.tvTran}
\title{Extract Coefficients from Time-varying Transformation Model}
\description{
  Extract coefficient values from \code{tvTran} fitting results, and
  summarize the point estimate and 95\% confidence band into a data
  frame.
}
\usage{
\method{coef}{tvTran}(object, \dots)
}
\arguments{
  \item{object}{
    an object returned by function \code{tvTran}.
  }
  \item{\dots}{
    other arguments.
  }
}
\value{
  A data.frame with 6 columns \code{("Low", "Mid", "High", "Time",
    "Cov", "Model")}, where \code{"Mid"} is the point estimates;
  \code{"Low"} and \code{"High"} are the 2.5\% and 97.5\% quantiles
  estimates from resampling method; \code{"Cov"} and \code{"Model"}
  contain character values of the covariates and model type.
}
\seealso{
  \code{\link{tvTran}}, and \code{\link{plotCoef}}.
}
\examples{
# See the examples in tvTran
}
\keyword{extract}
\keyword{tvTran}
\keyword{coefficient}

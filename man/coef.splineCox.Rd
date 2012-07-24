\name{coef.splineCox}
\alias{coef.splineCox}
\title{Extract Coefficients from Spline Base Cox Model}
\description{
  Extract coefficient values from \code{splineCox} fitting results, and
  summarize the point estimate and 95\% confidence band into a data
  frame.
}
\usage{
\method{coef}{splineCox}(object, \dots)
}
\arguments{
  \item{object}{
    an object returned by function \code{splineCox}.
  }
  \item{\dots}{
    other arguments.
  }
}
\value{
  A data.frame with 6 columns \code{("Low", "Mid", "High", "Time",
    "Cov", "Model")}, where \code{"Mid"} is the point estimates;
  \code{"Low"} and \code{"High"} are the point estimates plus and minus
  1.96 times standard deviations; \code{"Cov"} and \code{"Model"}
  contain character values of the covariates and model type.
}
\note{
  It essentially expand the break points, and then call function
  \code{coxph} in package \code{survival}
}
\seealso{
  \code{\link{splineCox}}, and \code{\link{plotCoef}}.
}
\examples{
# See the examples in splineCox
}
\keyword{extract}
\keyword{splineCox}
\keyword{coefficient}

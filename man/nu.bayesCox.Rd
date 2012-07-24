\name{nu.bayesCox}
\alias{nu}
\alias{nu.bayesCox}
\title{Extract Latent Variance from Bayesian Cox Model}
\description{
  Extract latent variance from \code{bayesCox} fitting results, and
  summarize them into a data frame. It is applicable when
  \code{model="TimeVarying"} or \code{model="Dynamic"}, and
  \code{coef.prior=list(type="HAR1")}.
}
\usage{
\method{nu}{bayesCox}(object, \dots)
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
  A data.frame with 4 columns \code{("Iter", "Model", "Cov", "Value")},
  where \code{Iter} is the iteration number; \code{Model} and \code{Cov}
  contain the character values of the model type and covariates.
}
\seealso{
  \code{\link{bayesCox}}, and \code{\link{plotNu}}.
}
\examples{
# See the examples in bayesCox
}
\keyword{extract}
\keyword{bayesCox}
\keyword{latent variance}
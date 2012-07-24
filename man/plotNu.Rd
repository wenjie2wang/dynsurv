\name{plotNu}
\alias{plotNu}
\title{Plot Latent Variance in Bayesian Cox Model}
\description{
  Plot the latent variance \code{nu} when the hierarchical AR(1) process
  prior is used for the \code{bayesCox} model. It is applicable when
  \code{model="TimeVarying"} or \code{model="Dynamic"}, and
  \code{coef.prior=list(type="HAR1")}. The input data frame is returned
  by function \code{nu}.
}
\usage{
plotNu(object, \dots)
}
\arguments{
  \item{object}{
    a data.frame returned by the function \code{nu}.
  }
  \item{\dots}{
    other arguments.
  }
}
\value{
  A \code{ggplot} object.
}
\seealso{
  \code{\link{nu.bayesCox}}.
}
\examples{
# See the examples in bayesCox
}
\keyword{plot}
\keyword{latent variance}

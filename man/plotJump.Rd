\name{plotJump}
\alias{plotJumpTrace}
\alias{plotJumpHist}
\title{Plot Jump Information in Bayesian Dynamic Model}
\description{
  \code{plotJumpTrace} plots the MCMC history of the number of
  pieces. \code{plotJumpHist} plots the histogram of the number of
  pieces. The input data frame is returned by function \code{jump}.
}
\usage{
plotJumpTrace(object, \dots)
plotJumpHist(object, \dots)
}
\arguments{
  \item{object}{
    a data.frame returned by function \code{jump}.
  }
  \item{\dots}{
    other arguments.
  }
}
\value{
  A \code{ggplot} object.
}
\seealso{
  \code{\link{jump.bayesCox}}.
}
\examples{
# See the examples in bayesCox
}
\keyword{plot}
\keyword{jump}

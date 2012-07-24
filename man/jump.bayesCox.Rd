\name{jump.bayesCox}
\alias{jump}
\alias{jump.bayesCox}
\title{Extract Jump Information from Bayesian Dyanmic Model}
\description{
  Extract number of coefficient pieces from \code{bayesCox} fitting
  results, and summarize them into a data frame. It is only applicable
  when \code{model="Dynamic"} is specified.
}
\usage{
\method{jump}{bayesCox}(object, \dots)
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
  A data.frame with 3 columns \code{("Count", "Iter", "Cov")}, where
  \code{"Count"} is the number of coefficient pieces (jumps) for each
  iteration; \code{Iter} is the iteration number; \code{Cov} contains
  the character values of the covariates.
}
\seealso{
  \code{\link{bayesCox}}, \code{\link{plotJumpTrace}}, and
  \code{\link{plotJumpHist}}.
}
\examples{
# See the examples in bayesCox
}
\keyword{extract}
\keyword{bayesCox}
\keyword{jump}

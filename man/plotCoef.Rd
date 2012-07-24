\name{plotCoef}
\alias{plotCoef}
\title{Plot Coefficient Function}
\description{
  Plot coefficient values formatted in a data frame returned by
  function \code{coef}.
}
\usage{
plotCoef(object, smooth=FALSE, \dots)
}
\arguments{
  \item{object}{
    a data.frame returned by function \code{coef}.
  }
  \item{smooth}{
    a logical value, default \code{FALSE}. If \code{TRUE}, plot the
    coefficients as smooth lines; otherwise, plot the coefficients as
    piece-wise constant step functions.
  }
  \item{\dots}{
    other arguments.
  }
}
\details{
  To plot estimated coefficient functions from different models
  together, one can first combine the data frames returned by
  \code{coef}, and then call \code{plotCoef}, for example,
  \code{plotCoef(rbind(coef(fit1), coef(fit2)))}.

  To specify the time range of the plot, one can either utilize the
  \code{ggplot} functionality, say \code{plotCoef(coef(fit)) +
    xlim(2, 10)}; or manipulate the data frame first before calling
  \code{plotCoef}, e.g., \code{plotCoef(subset(coef(fit), Time > 2 &
  Time < 10))}.
}
\value{
  A \code{ggplot} object.
}
\seealso{
  \code{\link{coef.bayesCox}}, \code{\link{coef.splineCox}}, and \code{\link{coef.tvTran}}.
}
\examples{
# See the examples in bayesCox, splineCox, and tvTran.
}
\keyword{plot}
\keyword{coefficient}

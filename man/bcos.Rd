\name{bcos}
\alias{bcos}
\alias{bcos.grid}
\docType{data}
\title{Breast Cancer Data}
\description{
  The breast cancer data from Finkelstein (1985) has been analyzed
  extensively for illustrating new methods in modeling interval
  censored data. The objective of the study was to compare the time to
  cosmetic deterioration between two groups: 46 patients receiving
  radiotherapy only and 48 patients receiving radiotherapy plus
  chemotherapy. Because the detection of deterioration required a
  clinic visit, the 56 women who experience deterioration were interval
  censored, and the other38 women who did not were right censored.
}
\usage{
bcos
bcos.grid
}
\format{
  \code{bcos} is a data frame with 94 observations and 3 columns.
  \describe{
    \item{\code{left}:}{left censoring time.}
    \item{\code{right}:}{right censoring time.}
    \item{\code{trt}:}{treatment, (\code{Rad} = radiotherapy only,
      \code{RadChem} = radiotherapy plus chemotherapy).}
  }

  \code{bcos.grid} is a numeric vector of grid time points.
}
%\details{}
%\source{
%}
\references{
  D.M. Finkelstein, and R.A. Wolfe (1985). A semiparametric model for
  regression analysis of interval-censored failure time
  data. \emph{Biometrics} 41: 731-740.
}
\examples{
data(bcos)
}
\keyword{datasets}

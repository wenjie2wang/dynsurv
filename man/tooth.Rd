\name{tooth}
\alias{tooth}
\alias{tooth.grid}
\docType{data}
\title{Tooth Data}
\description{
  The tooth data was from a longitudinal prospective dental study
  performed in Flanders (Belgium) in 1996 -- 2001. Every one of
  4,386 randomly sampled children in the cohort was examined annually
  by one of 16 trained dentists, resulting at most 6 dental observations
  for each child. The outcome of interest was the time to emergence of
  permanent tooth 24, which was either interval censored (2,775, 63\%)
  or right censored (1,611, 37\%).
}
\usage{
tooth
tooth.grid
}
\format{
  \code{tooth} is a data frame with 4,386 observations and 7 columns
  \describe{
    \item{\code{id}:}{children's id.}
    \item{\code{left}:}{left censoring time.}
    \item{\code{right}:}{right censoring time where infinity is coded as
      999.}
    \item{\code{sex}:}{gender of children (0 = boy, 1 = girl).}
    \item{\code{dmf}:}{status of the primary predecessor of this tooth
      (0 = sound, 1 = delayed, missing or filled).}
    \item{\code{rightInf}:}{right censoring time where infinity is coded
      as \code{Inf}.}
    \item{\code{rightNA}:}{right censoring time where infinity is coded
      as \code{NA}.}
  }

  \code{tooth.grid} is a numeric vector of grid time points.
}
%\details{}
\source{
  Adapted from the data set available at
  \url{http://www-eio.upc.es/grass/tips/8_tooth24.RData}.
}
\references{
  J. Vanobbergen, L. Martens, D. Declerck, and M. Lesaffre (2000). The
  signal tandmobiel(r) project: a longitudinal intervention oral health
  promotion study in Flanders (Belgium): baseline and first year
  results. \emph{European Journal of Paediatric Dentistry} 2, 87.

  G. Gomez, M. Calle, R. Oller, and K. Langohr (2009). Tutorial on
  methods for interval-censored data and their implementation in
  R. \emph{Statistical Modeling} 9(4), 259.
}
\examples{
data(tooth)
}
\keyword{datasets}

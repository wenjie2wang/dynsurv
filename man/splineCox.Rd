\name{splineCox}
\alias{splineCox}
\title{Fit Spline Based Cox Model for Right Censored Survival Data}
\description{
  Rearrange the rignt censored survival data in a counting process
  style. Model the time-varying coefficient function using
  B-splines. The fit is done by introducing pseudo time-dependent
  covariates and then calling function \code{coxph} in \code{survival}
  package.
}
\usage{
splineCox(formula, data, control=list())
}
\arguments{
  \item{formula}{
    a formula object, with the response on the left of a
    '~' operator, and the terms on the right. The response must be a
    survival object as returned by the \code{Surv} function.
  }
  \item{data}{
    a data.frame in which to interpret the variables named in
    the \code{formula}.
  }
  \item{control}{
    list of control options.
  }
}
\details{
  The \code{control} argument is a list of components:
  \describe{
    \item{\code{df}:}{degree of freedom for the B-splines, default 5.}
    \item{\code{knots}:}{interior knots point, default \code{NULL}. If
      \code{NULL}, the knots will be automatically choosen.}
    \item{\code{boundary}:}{lower and upper boundaries for the spline
      function, default \code{NULL}. If \code{NULL}, the minimun and maximun
      finite event time or censoring time will be specified.}
  }
}
\value{
  An object of S3 class \code{splineCox} representing the fit.
}
\references{
  A. Perperoglou, S. le Cessie, and H.C. van Houwelingen (2006). A fast
  routine for fitting Cox models with time varying effects of the
  covariates. \emph{Computer Methods and Programs in Biomedicine} 81,
  154--161.
}
\note{
  This function is essentially a wrapper function of \code{coxph} for the
  expanded data set. It does not implements the algorithm disscussed in the
  reference paper. These authors implemented their algorithm into a
  \code{tvcox} package, which is more efficient for larger data set, but
  may not be stable compared to \code{coxph}.
}
\seealso{
  \code{\link{coef.splineCox}}, \code{\link{plotCoef}}.
}
\examples{
\dontrun{
# Load the veteran data from the survival package
data <- veteran
data$celltype <- relevel(data$celltype, ref="large")
formula <- Surv(time, status) ~ karno + celltype

# Fit the time-varying transformation model
fit <- splineCox(formula, data, control=list(df=5))

# Plot the time-varying coefficient function between two time points
plotCoef(subset(coef(fit), Time > 15 & Time < 175), smooth=TRUE)
}
}
\keyword{B-spline}
\keyword{Cox}
\keyword{right censor}

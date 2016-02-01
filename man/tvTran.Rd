\name{tvTran}
\alias{tvTran}
\title{Fit Time-varying Transformation Model for Right Censored Survival
  Data}
\description{
  Unlike the time-varying coefficient Cox model, the transformation
  model fomulates the temporal covariate effects in terms of survival
  function, i.e., \deqn{S(t|X) = g(\beta_0(t)' X),} where \eqn{g(z) =
    exp(- exp(z))}. It can be viewed as a functional generalized linear
  model with response \eqn{I(T > t)}, and other transformation function
  is possible. The time-varying coefficients are solved a set of
  estimating equations sequentially.
}
\usage{
tvTran(formula, data, control=list())
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
  Note that because the time-varying coefficient function is connected
  to the survival function, it has a different interpretation of the
  time-varying coefficient function in Cox model.

  The \code{control} argument is a list of components:
  \describe{
    \item{\code{resample}:}{a logical value, default \code{TRUE}. If
      \code{TRUE}, the model will estimate a 95\% confidence band by
      resampling method.}
    \item{\code{R}:}{number of resamplings, default 30.}
  }
}
\value{
  An object of S3 class \code{tvTran} representing the fit.
}
\references{
  L. Peng, and Y. Huang (2007). Survival analysis with temporal
  covariate effects. \emph{Biometrika} 94(3), 719--733.
}
\seealso{
  \code{\link{coef.tvTran}}, \code{\link{plotCoef}}.
}
\examples{
\dontrun{
# Load the veteran data from the survival package
mydata <- survival::veteran
mydata$celltype <- relevel(mydata$celltype, ref="large")
myformula <- Surv(time, status) ~ karno + celltype

# Fit the time-varying transformation model
fit <- tvTran(myformula, mydata, control=list(resample=TRUE, R=30))

# Plot the time-varying coefficient function between two time points
plotCoef(subset(coef(fit), Time > 15 & Time < 175))
}
}
\keyword{transformation}
\keyword{right censor}

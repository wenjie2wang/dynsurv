\name{bayesCox}
\alias{bayesCox}
\title{Fit Bayesian Cox Model for Interval Censored Survival Data}
\description{
  Fit Bayesian Cox model with time-independent, time-varying or dynamic
  covariate coefficient. The fit is done within a Gibbs sampling
  framework. The reversible jump algorithm is employed for the dynamic
  coefficient model. The baseline hazards are allowed to be either
  time-varying or dynamic.
}
\usage{
bayesCox(formula, data, grid, out,
         model=c("TimeIndep", "TimeVarying", "Dynamic"),
         base.prior=list(), coef.prior=list(),
         gibbs=list(), control=list())
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
  \item{grid}{
    vector of pre-specified time grid points.
  }
  \item{out}{
    name of Markov chain Monte Carlo (MCMC) samples output file.
  }
  \item{model}{
    model type to fit.
  }
  \item{base.prior}{
    list of options for prior of baseline lambda. Use
    \code{list(type="Gamma", shape=0.1, rate=0.1)} for all models;
    \code{list(type="Const", value=1)} for \code{Dynamic} model when
    \code{intercept=TRUE}.
  }
  \item{coef.prior}{
    list of options for prior of coefficient beta. Use
    \code{list(type="Normal", mean=0, sd=1)} for \code{TimeIndep}
    model; \code{list(type="AR1", sd=1)} for \code{TimeVarying} and
    \code{Dynamic} models; \code{list(type="HAR1", shape=2, scale=1)} for
    \code{TimeVarying} and \code{Dynamic} models.
  }
  \item{gibbs}{
    list of options for Gibbs sampler.
  }
  \item{control}{
    list of general control options.
  }
}
\details{
  To use default hyper parameters in the specification of either
  \code{base.prior} or \code{coef.prior}, one only has to supply the
  name of the prior, e.g., \code{list(type="Gamma")},
  \code{list(type="HAR1")}.

  The \code{gibbs} argument is a list of components:
  \describe{
    \item{\code{iter}:}{number of iterations, default 3000.}
    \item{\code{burn}:}{number of burning, default 500.}
    \item{\code{thin}:}{number of thinning, default 1.}
    \item{\code{verbose}:}{a logical value, default \code{TRUE}. If
      \code{TRUE}, print the iteration.}
    \item{\code{nReport}:}{print frequency, default 100.}
  }

  The \code{control} argument is a list of components:
  \describe{
    \item{\code{intercept}:}{a logical value, default \code{FALSE}. If
      \code{TRUE}, the model will estimate the intercept, which is the
      log of baseline hazards. If \code{TRUE}, please remember to turn
      off the direct estimation of baseline hazards, i.e.,
      \code{base.prior=list(type="Const")}.}
    \item{\code{a0}:}{multiplier for initial variance in time-varying or
      dynamic models, default 100.}
    \item{\code{eps0}:}{size of auxiliary uniform latent variable in
      dynamic model, default 1.}
  }
}
\value{
  An object of S3 class \code{bayesCox} representing the fit.
}
\references{
  X. Wang, M.-H. Chen, and J. Yan (2011). Bayesian dynamic regression
  models for interval censored survival data. Under review.

  D. Sinha, M.-H. Chen, and S.K. Ghosh (1999). Bayesian analysis and
  model selection for interval-censored survival data. \emph{Biometrics}
  55(2), 585--590.
}

\seealso{
  \code{\link{coef.bayesCox}}, \code{\link{jump.bayesCox}},
  \code{\link{nu.bayesCox}}, \code{\link{plotCoef}},
  \code{\link{plotJumpTrace}}, and \code{\link{plotNu}}.
}
\examples{
\dontrun{
################################################################################
# Load one of the following two data sets
################################################################################

# breast cancer data
data(bcos) ## load bcos and bcos.grid
mydata <- bcos
mygrid <- bcos.grid
myformula <- Surv(left, right, type="interval2") ~ trt

# tooth data
# data(tooth) ## load tooth and tooth.grid
# mydata <- tooth
# mygrid <- tooth.grid
# myformula <- Surv(left, right, type="interval2") ~ dmf + sex

################################################################################
# Fit Bayesian Cox models
################################################################################

# Fit time-independent coefficient model
fit0 <- bayesCox(myformula, mydata, mygrid, out="tiCox.txt",
                 model="TimeIndep",
                 base.prior=list(type="Gamma", shape=0.1, rate=0.1),
                 coef.prior=list(type="Normal", mean=0, sd=1),
                 gibbs=list(iter=100, burn=20, thin=1, verbose=TRUE, nReport=5))
plotCoef(coef(fit0))

# Fit time-varying coefficient model
fit1 <- bayesCox(myformula, mydata, mygrid, out="tvCox.txt",
                 model="TimeVarying",
                 base.prior=list(type="Gamma", shape=0.1, rate=0.1),
                 coef.prior=list(type="AR1", sd=1),
                 gibbs=list(iter=100, burn=20, thin=1, verbose=TRUE, nReport=5))
plotCoef(coef(fit1))

# Fit dynamic coefficient model with time-varying baseline hazards
fit2 <- bayesCox(myformula, mydata, mygrid, out="dynCox1.txt",
                 model="Dynamic",
                 base.prior=list(type="Gamma", shape=0.1, rate=0.1),
                 coef.prior=list(type="HAR1", shape=2, scale=1),
                 gibbs=list(iter=100, burn=20, thin=1, verbose=TRUE, nReport=5))
plotCoef(coef(fit2))
plotJumpTrace(jump(fit2))
plotJumpHist(jump(fit2))
plotNu(nu(fit2))

# Plot the coefficient estimates from three models together
plotCoef(rbind(coef(fit0), coef(fit1), coef(fit2)))

# Fit dynamic coefficient model with dynamic hazards (in log scales)
fit3 <- bayesCox(myformula, mydata, mygrid, out="dynCox2.txt",
                 model="Dynamic",
                 base.prior=list(type="Const"),
                 coef.prior=list(type="HAR1", shape=2, scale=1),
                 gibbs=list(iter=100, burn=20, thin=1, verbose=TRUE, nReport=5),
                 control=list(intercept=TRUE))
plotCoef(coef(fit3))
plotJumpTrace(jump(fit3))
plotJumpHist(jump(fit3))
plotNu(nu(fit3))
}
}
\keyword{Bayesian}
\keyword{Cox}
\keyword{dynamic}
\keyword{interval censor}

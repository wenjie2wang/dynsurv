################################################################################
##
##   R package dynsurv by Wenjie Wang, Ming-Hui Chen, Xiaojing Wang, and Jun Yan
##   Copyright (C) 2011-2016
##
##   This file is part of the R package dynsurv.
##
##   The R package dynsurv is free software: You can redistribute it and/or
##   modify it under the terms of the GNU General Public License as published
##   by the Free Software Foundation, either version 3 of the License, or
##   any later version (at your option). See the GNU General Public License
##   at <http://www.gnu.org/licenses/> for details.
##
##   The R package dynsurv is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##
################################################################################


### Bayesian Cox model
## grid: must be sorted with last number be finite
## base.prior:
##   list(type = "Gamma", shape = 0.1, rate = 0.1)
## coef.prior:
##   list(type = "Normal", mean = 0, sd = 1)
##   list(type = "AR1", sd = 1)
##   list(type = "HAR1", shape = 2, scale = 1)

##' Fit Bayesian Cox Model for Interval Censored Survival Data
##'
##' Fit Bayesian Cox model with time-independent, time-varying or dynamic
##' covariate coefficient. The fit is done within a Gibbs sampling framework.
##' The reversible jump algorithm is employed for the dynamic coefficient
##' model. The baseline hazards are allowed to be either time-varying or
##' dynamic.
##'
##' For application, it is recommended that users preprocess the data by
##' rounding down (up) the left (right) endpoints of the censoring intervals
##' to at most three significant digits and leave the \code{grid} unspecified
##' in the function call to reduce probably unnecessary computational burden.
##' It also helps avoid possible errors caused by mispecifed \code{grid}. For
##' example, the left endpoints can be rounded down in the nearest 0.1 unit
##' of time by \code{left <- floor(left * 10) / 10} and the right endpoints
##' can similarly be rounded up by \code{right <- ceiling(right * 10) / 10}.
##'
##' To use default hyper parameters in the specification of either
##' \code{base.prior} or \code{coef.prior}, one only has to supply the name of
##' the prior, e.g., \code{list(type = "Gamma")}, \code{list(type = "HAR1")}.
##'
##' The \code{gibbs} argument is a list of components:
##' \describe{
##'     \item{iter:}{number of iterations, default 3000;}
##'     \item{burn:}{number of burning, default 500;}
##'     \item{thin:}{number of thinning, default 1;}
##'     \item{verbose:}{a logical value, default \code{TRUE}. If
##'         \code{TRUE}, print the iteration;}
##'     \item{nReport:}{print frequency, default 100.}
##' }
##'
##' The \code{control} argument is a list of components:
##' \describe{
##'     \item{intercept:}{a logical value, default \code{FALSE}. If
##'         \code{TRUE}, the model will estimate the intercept, which is the
##'         log of baseline hazards. If \code{TRUE}, please remember to turn
##'         off the direct estimation of baseline hazards, i.e.,
##'         \code{base.prior = list(type = "Const")}}
##'     \item{a0:}{multiplier for initial variance in time-varying or dynamic
##'         models, default 100;}
##'     \item{eps0:}{size of auxiliary uniform latent variable in dynamic model,
##'         default 1.}
##' }
##'
##' @usage
##' bayesCox(formula, data, grid = NULL, out = "mcmc.txt", model =
##'          c("TimeIndep", "TimeVarying", "Dynamic"), base.prior = list(),
##'          coef.prior = list(), gibbs = list(), control = list())
##'
##' @param formula A formula object, with the response on the left of a '~'
##'     operator, and the terms on the right. The response must be a survival
##'     object as returned by the \code{Surv} function.
##' @param data A data.frame in which to interpret the variables named in the
##'     \code{formula}.
##' @param grid Vector of pre-specified time grid points for model fitting.
##'     It will be automatically set up from data if it is left unspecified
##'     in the function call. By default, it consists of all the unique
##'     finite endpoints of the censoring intervals after time zero.  The
##'     \code{grid} specified in the function call must be sorted, and covers
##'     all the finite non-zero endpoints of the censoring intervals.
##' @param out Name of Markov chain Monte Carlo (MCMC) samples output file.
##'     Each row contains one MCMC sample information. The file is needed for
##'     those functions further summarizing estimation results in this
##'     package.
##' @param model Model type to fit.
##' @param base.prior List of options for prior of baseline lambda. Use
##'     \code{list(type = "Gamma", shape = 0.1, rate = 0.1)} for all models;
##'     \code{list(type = "Const", value = 1)} for \code{Dynamic} model when
##'     \code{intercept = TRUE}.
##' @param coef.prior List of options for prior of coefficient beta. Use
##'     \code{list(type = "Normal", mean = 0, sd = 1)} for \code{TimeIndep}
##'     model; \code{list(type = "AR1", sd = 1)} for \code{TimeVarying} and
##'     \code{Dynamic} models; \code{list(type = "HAR1", shape = 2, scale =
##'     1)} for \code{TimeVarying} and \code{Dynamic} models.
##' @param gibbs List of options for Gibbs sampler.
##' @param control List of general control options.
##' @return An object of S3 class \code{bayesCox} representing the fit.
##' @seealso \code{\link{coef.bayesCox}}, \code{\link{jump.bayesCox}},
##'     \code{\link{nu.bayesCox}}, \code{\link{plotCoef}},
##'     \code{\link{plotJumpTrace}}, \code{\link{plotNu}},
##'     \code{\link{survCurve}}, \code{\link{survDiff}}, and
##'     \code{\link{plotSurv}}.
##'
##' @references
##' X. Wang, M.-H. Chen, and J. Yan (2013). Bayesian dynamic regression
##' models for interval censored survival data with application to children
##' dental health. Lifetime data analysis, 19(3), 297--316.
##'
##' X. Wang, X. Sinha, J. Yan, and M.-H. Chen (2014). Bayesian inference of
##' interval-censored survival data. In: D. Chen, J. Sun, and K. Peace,
##' Interval-censored time-to-event data: Methods and applications, 167--195.
##'
##' X. Wang, M.-H. Chen, and J. Yan (2011). Bayesian dynamic
##' regression models for interval censored survival data. Technical Report 13,
##' Department of Statistics, University of Connecticut.
##'
##' D. Sinha, M.-H. Chen, and S.K. Ghosh (1999). Bayesian analysis and model
##' selection for interval-censored survival data. \emph{Biometrics} 55(2),
##' 585--590.
##' @keywords Bayesian Cox dynamic interval censor
##' @examples
##' \dontrun{
##' ############################################################################
##' ### Attach one of the following two data sets
##' ############################################################################
##'
##' ## breast cancer data
##' data(bcos) ## attach bcos and bcos.grid
##' mydata <- bcos
##' ## mygrid <- bcos.grid
##' myformula <- Surv(left, right, type = "interval2") ~ trt
##'
##' ## tooth data
##' ## data(tooth) ## load tooth and tooth.grid
##' ## mydata <- tooth
##' ## mygrid <- tooth.grid
##' ## myformula <- Surv(left, rightInf, type = "interval2") ~ dmf + sex
##'
##' ############################################################################
##' ### Fit Bayesian Cox models
##' ############################################################################
##'
##' ## Fit time-independent coefficient model
##' fit0 <- bayesCox(myformula, mydata, out = "tiCox.txt",
##'                  model = "TimeIndep",
##'                  base.prior = list(type = "Gamma", shape = 0.1, rate = 0.1),
##'                  coef.prior = list(type = "Normal", mean = 0, sd = 1),
##'                  gibbs = list(iter = 100, burn = 20, thin = 1,
##'                               verbose = TRUE, nReport = 5))
##' plotCoef(coef(fit0, level = 0.9))
##'
##' ## Fit time-varying coefficient model
##' fit1 <- bayesCox(myformula, mydata, out = "tvCox.txt",
##'                  model = "TimeVarying",
##'                  base.prior = list(type = "Gamma", shape = 0.1, rate = 0.1),
##'                  coef.prior = list(type = "AR1", sd = 1),
##'                  gibbs = list(iter = 100, burn = 20, thin = 1,
##'                               verbose = TRUE, nReport = 5))
##' plotCoef(coef(fit1))
##'
##' ## Fit dynamic coefficient model with time-varying baseline hazards
##' fit2 <- bayesCox(myformula, mydata, out = "dynCox1.txt",
##'                  model = "Dynamic",
##'                  base.prior = list(type = "Gamma", shape = 0.1, rate = 0.1),
##'                  coef.prior = list(type = "HAR1", shape = 2, scale = 1),
##'                  gibbs = list(iter = 100, burn = 20, thin = 1,
##'                               verbose = TRUE, nReport = 5))
##' plotCoef(coef(fit2))
##' plotJumpTrace(jump(fit2))
##' plotJumpHist(jump(fit2))
##' plotNu(nu(fit2))
##'
##' ## Plot the coefficient estimates from three models together
##' plotCoef(rbind(coef(fit0), coef(fit1), coef(fit2)))
##'
##' ## Fit dynamic coefficient model with dynamic hazards (in log scales)
##' fit3 <- bayesCox(myformula, mydata, out = "dynCox2.txt",
##'                  model = "Dynamic",
##'                  base.prior = list(type = "Const"),
##'                  coef.prior = list(type = "HAR1", shape = 2, scale = 1),
##'                  gibbs = list(iter = 100, burn = 20, thin = 1,
##'                               verbose = TRUE, nReport=5),
##'                  control = list(intercept = TRUE))
##' plotCoef(coef(fit3))
##' plotJumpTrace(jump(fit3))
##' plotJumpHist(jump(fit3))
##' plotNu(nu(fit3))
##'
##' ## Plot the estimated survival function and hazard function
##' newDat <- bcos[c(1L, 47L), ]
##' row.names(newDat) <- c("Rad", "RadChem")
##' plotSurv(survCurve(fit3, newdata = newDat, type = "survival"),
##'          legendName = "Treatment", conf.int = TRUE)
##' plotSurv(survDiff(fit3, newdata = newDat, type = "cumhaz"),
##'          legendName = "Treatment", conf.int = TRUE, smooth = TRUE)
##' }
##'
##' @importFrom stats model.frame model.matrix .getXlevels
##' @importFrom utils tail
##' @export bayesCox
bayesCox <- function(formula, data, grid = NULL, out = "mcmc.txt",
                     model = c("TimeIndep", "TimeVarying", "Dynamic"),
                     base.prior = list(), coef.prior = list(), gibbs = list(),
                     control = list()) {
    Call <- match.call()

    ## Prepare prior information
    model <- match.arg(model)
    base.prior <- do.call("bp_fun", base.prior)
    coef.prior <- do.call("cp_fun", coef.prior)

    if (model == "TimeIndep") {
        if (base.prior$type == "Gamma" && coef.prior$type == "Normal")
            id <- 11
        if (base.prior$type == "GammaProcess" &&
            coef.prior$type == "Normal")
            id <- 12
    }
    if (model == "TimeVarying") {
        if (base.prior$type == "Gamma" && coef.prior$type == "AR1")
            id <- 21
        if (base.prior$type == "Gamma" && coef.prior$type == "HAR1")
            id <- 22
        if (base.prior$type == "GammaProcess" &&
            coef.prior$type == "AR1")
            id <- 23
        if (base.prior$type == "GammaProcess" &&
            coef.prior$type == "HAR1")
            id <- 24

    }
    if (model == "Dynamic") {
        if (base.prior$type == "Gamma" && coef.prior$type == "AR1")
            id <- 31
        if (base.prior$type == "Gamma" && coef.prior$type == "HAR1")
            id <- 32
        if (base.prior$type == "Const" && coef.prior$type == "AR1")
            id <- 33
        if (base.prior$type == "Const" && coef.prior$type == "HAR1")
            id <- 34
        if (base.prior$type == "GammaProcess" &&
            coef.prior$type == "AR1")
            id <- 35
        if (base.prior$type == "GammaProcess" &&
            coef.prior$type == "HAR1")
            id <- 36
    }

    gibbs <- do.call("gibbs_fun", gibbs)
    control <- do.call("control_bfun", control)

    ## Prepare data matrix LRX
    mf <- model.frame(formula, data)
    mt <- attr(mf, "terms")
    mm <- model.matrix(formula, data)

    LRX <- cbind(mf[, 1][, 1:2], mm[, -1])
    obsInd <- which(mf[, 1][, 3] == 1)
    LRX[obsInd, 2] <- LRX[obsInd, 1]

    cov.names <- colnames(mm)[-1]

    ## prepare the grid if it is not specified
    if (is.null(grid)) {
        ## determine the significant digits from the data
        charNum <- as.character(LRX[, seq_len(2)])
        tmpList <- strsplit(charNum, "\\.")
        sigMax <- max(sapply(tmpList, function(a) {
            if (length(a) > 1)
                return(nchar(a[2L]))
            0
        }))
        if (sigMax > 3)
            warning(paste("Endpoints of censoring intervals were enlarged",
                          "to three significant digits\n to reduce",
                          "computational burden."))
        ## round the right (left) endpoints up (down)
        roundPow <- 10 ^ min(sigMax, 3)
        left_ <- floor(LRX[, "time1"] * roundPow) / roundPow
        right_ <- ceiling(LRX[, "time2"] * roundPow) / roundPow
        finiteRight <- right_[! (is.infinite(right_) | is.na(right_))]
        ## set up grid from the data
        grid <- sort(unique(c(left_, finiteRight)))
        ## make sure the grid does not contain zero
        grid <- grid[grid > 0]

        LRX[, "time1"] <- left_
        LRX[, "time2"] <- right_
    }

    LRX[LRX[, 2] == Inf, 2] <- max(tail(grid, 1), 999)
    LRX[is.na(LRX[, 2]), 2] <- max(tail(grid, 1), 999)

    if (control$intercept) {
        LRX <- cbind(LRX[, 1:2], 1, LRX[, -c(1:2)])
        cov.names <- c("intercept", cov.names)
    }

    colnames(LRX) <- c("L", "R", cov.names)

    ## Prepare results holder
    K <- length(grid)
    nBeta <- length(cov.names)
    lambda <- rep(0, K)

    if (model == "TimeIndep")
        beta <- rep(0, nBeta)
    else
        beta <- rep(0, nBeta * K)

    nu <- rep(0, nBeta)
    jump <- rep(0, nBeta * K)

    ## Call C++ function
    res <- .C("bayesCox",
              as.double(LRX),
              as.integer(nrow(LRX)),
              as.integer(ncol(LRX) - 2),
              as.double(grid),
              as.integer(length(grid)),
              as.character(out),
              as.integer(id),
              as.double(base.prior$hyper[[1]]),
              as.double(base.prior$hyper[[2]]),
              as.double(coef.prior$hyper[[1]]),
              as.double(coef.prior$hyper[[2]]),
              as.integer(gibbs$iter),
              as.integer(gibbs$burn),
              as.integer(gibbs$thin),
              as.integer(gibbs$verbose),
              as.integer(gibbs$nReport),
              as.double(control$a0),
              as.double(control$eps0),
              lambda = as.double(lambda),
              beta = as.double(beta),
              nu = as.double(nu),
              jump = as.integer(jump),
              LPML = as.double(0),
              DHat = as.double(0),
              DBar = as.double(0),
              pD = as.double(0),
              DIC = as.double(0))

    ## Post fit processing
    if (model != "TimeIndep")
        res$beta <- matrix(res$beta, K, nBeta)

    if (coef.prior$type != "HAR1")
        res$nu <- NULL

    if (model != "Dynamic")
        res$jump <- NULL
    else
        res$jump <- matrix(res$jump, K, nBeta)

    ## Return list
    rl <- list(call = Call, formula = formula, grid = grid, out = out,
               model = model, LRX = LRX, base.prior = base.prior,
               coef.prior = coef.prior, gibbs = gibbs,
               control = control, xlevels = .getXlevels(mt, mf),
               N = nrow(LRX), K = K, nBeta = nBeta, cov.names = cov.names,
               est = list(lambda = res$lambda, beta = res$beta, nu = res$nu,
                          jump = res$jump),
               measure = list(LPML = res$LPML, DHat = res$DHat, DBar = res$DBar,
                              pD = res$pD, DIC = res$DIC))
    class(rl) <- "bayesCox"
    rl
}


### Internal functions =========================================================

## Baseline prior
Gamma_fun <- function(shape = 0.1, rate = 0.1) {
    list(shape = shape, rate = rate)
}

Const_fun <- function(value = 1) {
    list(value = value, value = value)
}

GammaProcess_fun <- function(mean = 0.1, ctrl = 1) {
    list(mean = mean, ctrl = ctrl)
}

bp_fun <- function(type = c("Gamma", "Const", "GammaProcess"), ...) {
    type <- match.arg(type)

    if (type == "Gamma")
        hyper <- Gamma_fun(...)
    if (type == "Const")
        hyper <- Const_fun(...)
    if (type == "GammaProcess")
        hyper <- GammaProcess_fun(...)

    list(type = type, hyper = hyper)
}


## Coefficient prior
Normal_fun <- function(mean = 0, sd = 1) {
    list(mean = mean, sd = sd)
}
AR1_fun <- function(sd = 1) {
    list(sd = sd, sd = sd)
}
HAR1_fun <- function(shape = 2, scale = 1) {
    list(shape = shape, scale = scale)
}

cp_fun <- function(type = c("Normal", "AR1", "HAR1"), ...) {
    type <- match.arg(type)

    if (type == "Normal")
        hyper <- Normal_fun(...)
    if (type == "AR1")
        hyper <- AR1_fun(...)
    if (type == "HAR1")
        hyper <- HAR1_fun(...)

    list(type = type, hyper = hyper)
}


## Gibbs sampler control and general control
gibbs_fun <- function(iter = 3000, burn = 500, thin = 1, verbose = TRUE,
                      nReport = 100) {
    list(iter = iter,
         burn = burn,
         thin = thin,
         verbose = verbose,
         nReport = nReport)
}

control_bfun <- function(intercept = FALSE, a0 = 100, eps0 = 1) {
    list(intercept = intercept, a0 = a0, eps0 = eps0)
}

##
## R package dynsurv by Wenjie Wang, Ming-Hui Chen, Xiaojing Wang, and Jun Yan
## Copyright (C) 2011-2020
##
## This file is part of the R package dynsurv.
##
## The R package dynsurv is free software: You can redistribute it and/or
## modify it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package dynsurv is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##


### Bayesian Cox model
## grid: must be sorted with last number being finite
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
##' To use default hyper parameters in the specification of either
##' \code{base.prior} or \code{coef.prior}, one only has to supply the name of
##' the prior, e.g., \code{list(type = "Gamma")}, \code{list(type = "HAR1")}.
##'
##' The \code{gibbs} argument is a list of components:
##' \describe{
##'     \item{iter:}{Number of iterations, default 3000;}
##'     \item{burn:}{Number of burning, default 500;}
##'     \item{thin:}{Number of thinning, default 1;}
##'     \item{verbose:}{A logical value, default \code{TRUE}. If
##'         \code{TRUE}, print the iteration;}
##'     \item{nReport:}{Print frequency, default 100.}
##' }
##'
##' The \code{control} argument is a list of components:
##' \describe{
##'     \item{intercept:}{A logical value, default \code{FALSE}. If
##'         \code{TRUE}, the model will estimate the intercept, which is the
##'         log of baseline hazards. If \code{TRUE}, please remember to turn
##'         off the direct estimation of baseline hazards, i.e.,
##'         \code{base.prior = list(type = "Const")}}
##'     \item{a0:}{Multiplier for initial variance in time-varying or dynamic
##'         models, default 100;}
##'     \item{eps0:}{Size of auxiliary uniform latent variable in dynamic model,
##'         default 1.}
##' }
##'
##' For users interested in extracting MCMC sampling information from the
##' output files, the detail of the output files is presented as follows: Let
##' \eqn{k} denote the number of time points (excluding time zero) specified
##' in grid, \eqn{ck} equal \eqn{1} for model with time-invariant coefficients;
##' \eqn{ck} equal \eqn{k} otherwise, and \eqn{p} denote the number of
##' covariates.  Then the each sample saved in each row consists of the
##' following possible parts.
##' \describe{
##'     \item{Part 1:}{The first \eqn{k} numbers represent the jump size of
##' baseline hazard function at each time grid.  If we take the column mean
##' of the first \eqn{k} columns of the output file, we will get the same
##' numbers with \code{obj$est$lambda}, where \code{obj} is the \code{bayesCox}
##' object returned by the function.}
##'     \item{Part 2:}{The sequence from \eqn{(k + 1) to (k + ck * p)}
##' represent the coefficients of covariates at the time grid.  The first
##' \eqn{k} numbers in the sequence are the coefficients for the first covariate
##' at the time grid; The second \eqn{k} numbers' sub-sequence are the
##' coefficients for the second covariate and so on.}
##'     \item{Part 3:}{The sequence from \eqn{(k + ck * p + 1)} to
##' \eqn{(k + ck * p + p)} represents the sampled latent variance of
##' coefficients.}
##'     \item{Part 4:}{The sequence from \eqn{(k + ck * p + p + 1)} to
##' \eqn{(k + 2 * ck * p + p)} represents the indicator of whether there is
##' a jump of the covariate coefficients at the time grid.  Similar with Part 2,
##' the first k numbers' sub-sequence is for the first covariate, the second
##' \eqn{k} numbers' sub-sequence is for the second covariate, and so on.}
##' }
##' For the model with time-independent coefficients, the output file only
##' has Part 1 and Part 2 in each row; For time-varying coefficient model,
##' the output file has Part 1, 2, and 3; The output file for the dynamic
##' model has all the four parts.  Note that the dynamic baseline hazard will
##' be taken as one covariate.  So \eqn{p} needs being replaced with
##' \eqn{(p + 1)} for model with dynamic baseline hazard rate.
##' No function in the package actually needs the Part 1 from the output file
##' now; The Part 2 is used by function \code{coef} and \code{survCurve};
##' The Part 3 is needed by function \code{nu}; Function \code{jump} extracts
##' the Part 4.
##'
##' @param formula A formula object, with the response on the left of a '~'
##'     operator, and the terms on the right. The response must be a survival
##'     object as returned by the function \code{Surv} with \code{type =
##'     "interval2"}. \code{help(Surv)} for details.
##' @param data A data.frame in which to interpret the variables named in the
##'     \code{formula}.
##' @param grid Vector of pre-specified time grid points for model fitting.  It
##'     will be automatically set up from data if it is left unspecified in the
##'     function call. By default, it consists of all the unique finite
##'     endpoints (rounded to two significant numbers) of the censoring
##'     intervals after time zero.  The \code{grid} specified in the function
##'     call determines the location of possible jumps. It should be sorted
##'     increasingly and cover all the finite non-zero endpoints of the
##'     censoring intervals. Inappropriate \code{grid} specified will be taken
##'     care by the function internally.
##' @param out An optional character string specifying the name of Markov chain
##'     Monte Carlo (MCMC) samples output file.  By default, the MCMC samples
##'     will be output to a temporary directory set by \code{tempdir} and saved
##'     in the returned \code{bayesCox} object after burning and thinning.  If
##'     \code{out} is specified, the MCMC samples will be preserved in the
##'     specified text file.
##' @param model Model type to fit. Available options are \code{"TimeIndep"},
##'     \code{"TimeVarying"}, and \code{"Dynamic"}. Partial matching on the name
##'     is allowed.
##' @param base.prior List of options for prior of baseline lambda. Use
##'     \code{list(type = "Gamma", shape = 0.1, rate = 0.1)} for all models;
##'     \code{list(type = "Const", value = 1)} for \code{Dynamic} model when
##'     \code{intercept = TRUE}.
##' @param coef.prior List of options for prior of coefficient beta. Use
##'     \code{list(type = "Normal", mean = 0, sd = 1)} for \code{TimeIndep}
##'     model; \code{list(type = "AR1", sd = 1)} for \code{TimeVarying} and
##'     \code{Dynamic} models; \code{list(type = "HAR1", shape = 2, scale = 1)}
##'     for \code{TimeVarying} and \code{Dynamic} models.
##' @param gibbs List of options for Gibbs sampler.
##' @param control List of general control options.
##' @param ... Other arguments that are for futher extension.
##'
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
##'
##' @keywords Bayesian Cox dynamic interval censor
##'
##' @example inst/examples/ex-bayesCox.R
##'
##' @importFrom stats model.frame model.matrix .getXlevels stepfun
##'
##' @export
bayesCox <- function(formula, data, grid = NULL, out = NULL,
                     model = c("TimeIndep", "TimeVarying", "Dynamic"),
                     base.prior = list(),
                     coef.prior = list(),
                     gibbs = list(),
                     control = list(),
                     ...)
{
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
    } else if (model == "TimeVarying") {
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
    } else if (model == "Dynamic") {
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
    } else
        stop("Invalid 'model' specified.")

    gibbs <- do.call("gibbs_fun", gibbs)
    control <- do.call("control_bfun", control)

    ## Prepare data matrix LRX
    mf <- model.frame(formula, data)
    mt <- attr(mf, "terms")
    mm <- model.matrix(formula, data)

    LRX <- cbind(mf[, 1L][, seq_len(2L)], mm[, - 1L])

    ## reference: Surv.S from survival package
    eventIdx <- mf[, 1L][, 3L]
    ## for possibly - Inf or NA left
    tmpIdx <- eventIdx == 2
    LRX[tmpIdx, 1L] <- 0

    ## for right censoring
    tmpIdx <- eventIdx == 0
    LRX[tmpIdx, 2L] <- Inf

    ## for exact event times
    tmpIdx <- eventIdx == 1
    exactTimes <- LRX[tmpIdx, 2L] <- LRX[tmpIdx, 1L]

    ## record coveriate names
    cov.names <- colnames(mm)[- 1L]

    ## take care of grid
    if (is.null(grid)) {
        ## prepare the grid if it is not specified
        ## determine the significant digits from the data
        charNum <- as.character(LRX[, seq_len(2)])
        tmpList <- strsplit(charNum, "\\.")
        sigMax <- max(sapply(tmpList, function(a) {
            if (length(a) > 1)
                return(nchar(a[2L]))
            0
        }))
        ## round the right (left) endpoints up (down)
        roundPow <- 10 ^ min(sigMax, 2)
        left_ <- floor(LRX[, "time1"] * roundPow) / roundPow
        right_ <- ceiling(LRX[, "time2"] * roundPow) / roundPow
        finiteRight <- right_[! (is.infinite(right_) | is.na(right_))]
        ## set up grid from the data
        grid <- unique(c(left_, finiteRight))
    }

    ## prepare data based on the grid
    ## make sure all grid points great than zero, finite and sorted
    grid <- grid[! (is.na(grid) | is.infinite(grid))]
    ## exclude non-positive points
    grid <- grid[! grid <= 0]
    ## add exact event times to the grid
    if (length(exactTimes) > 0) {
        grid <- unique(c(grid, exactTimes))
    }
    if (all(tmpIdx <- is.infinite(LRX[, "time2"]))) {
        stop("Subjects are all right censored.")
    }
    finiteRight <- max(LRX[! tmpIdx, "time2"])
    if (max(grid) < finiteRight) {
        warning("The grid was expanded to cover all the finite endpoint ",
                "of censoring intervals.")
        grid <- unique(c(grid, finiteRight))
    }

    ## make sure grid is sorted and consists of unique points
    grid <- sort(unique(grid))

    ## round the left (right) endpoints down (up)
    ## to the closest grid point including zero
    toLeft <- stats::stepfun(grid, c(0, grid))
    toRight <- stats::stepfun(grid, c(grid, Inf))
    LRX[, "time1"] <- toLeft(LRX[, "time1"])
    LRX[, "time2"] <- toRight(LRX[, "time2"])

    LRX[is.infinite(LRX[, 2L]), 2L] <- max(grid[length(grid)], 999)
    LRX[is.na(LRX[, 2L]), 2L] <- max(grid[length(grid)], 999)

    if (control$intercept) {
        LRX <- cbind(LRX[, seq_len(2L)], 1, LRX[, - seq_len(2L)])
        cov.names <- c("intercept", cov.names)
    }
    colnames(LRX) <- c("L", "R", cov.names)

    ## if out is not specified, use the temp file
    if (save_mcmc <- is.null(out)) {
        out <- sprintf("%s.txt", tempfile())
    }

    ## Prepare results holder
    K <- length(grid)
    nBeta <- length(cov.names)
    lambda <- rep(0, K)
    beta <- if (model == "TimeIndep") {
                rep(0, nBeta)
            } else {
                rep(0, nBeta * K)
            }
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

    if (model == "Dynamic") {
        res$jump <- matrix(res$jump, K, nBeta)
    } else {
        res$jump <- NULL
    }

    ## save the MCMC samples if out was not specified
    if (save_mcmc) {
        mcmc_list <- read_bayesCox(out, gibbs$burn, gibbs$thin)
        out <- NULL
    } else {
        mcmc_list <- NULL
    }

    ## Return a list of class bayesCox
    structure(list(
        call = Call,
        formula = formula,
        grid = grid,
        out = out,
        model = model,
        LRX = LRX,
        base.prior = base.prior,
        coef.prior = coef.prior,
        gibbs = gibbs,
        control = control,
        xlevels = .getXlevels(mt, mf),
        N = nrow(LRX),
        K = K,
        nBeta = nBeta,
        cov.names = cov.names,
        mcmc = mcmc_list,
        est = list(lambda = res$lambda,
                   beta = res$beta,
                   nu = res$nu,
                   jump = res$jump),
        measure = list(LPML = res$LPML,
                       DHat = res$DHat,
                       DBar = res$DBar,
                       pD = res$pD,
                       DIC = res$DIC)
    ), class = "bayesCox")
}


### Internal functions =========================================================

## Baseline prior
Gamma_fun <- function(shape = 0.1, rate = 0.1)
{
    list(shape = shape, rate = rate)
}

Const_fun <- function(value = 1)
{
    list(value = value, value = value)
}

GammaProcess_fun <- function(mean = 0.1, ctrl = 1)
{
    list(mean = mean, ctrl = ctrl)
}

bp_fun <- function(type = c("Gamma", "Const", "GammaProcess"), ...)
{
    type <- match.arg(type)
    hyper <- if (type == "Gamma") {
                 Gamma_fun(...)
             } else if (type == "Const") {
                 Const_fun(...)
             } else if (type == "GammaProcess") {
                 GammaProcess_fun(...)
             }
    list(type = type, hyper = hyper)
}


## Coefficient prior
Normal_fun <- function(mean = 0, sd = 1)
{
    list(mean = mean, sd = sd)
}
AR1_fun <- function(sd = 1)
{
    list(sd = sd, sd = sd)
}
HAR1_fun <- function(shape = 2, scale = 1)
{
    list(shape = shape, scale = scale)
}

cp_fun <- function(type = c("Normal", "AR1", "HAR1"), ...)
{
    type <- match.arg(type)
    hyper <- if (type == "Normal") {
                 Normal_fun(...)
             } else if (type == "AR1") {
                 AR1_fun(...)
             } else if (type == "HAR1") {
                 HAR1_fun(...)
             }
    list(type = type, hyper = hyper)
}


## Gibbs sampler control and general control
gibbs_fun <- function(iter = 3000, burn = 500, thin = 1,
                      verbose = TRUE, nReport = 100)
{
    list(iter = as.integer(iter),
         burn = as.integer(burn),
         thin = as.integer(thin),
         verbose = verbose,
         nReport = as.integer(nReport))
}

control_bfun <- function(intercept = FALSE, a0 = 100, eps0 = 1)
{
    list(intercept = intercept, a0 = a0, eps0 = eps0)
}

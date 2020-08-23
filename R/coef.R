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


##' Extract Coefficients from Bayesian Cox Model
##'
##' Extract coefficient values from \code{bayesCox} fitting results, and
##' summarize the posterior mean, posterior 2.5\% and 97.5\% quantiles into a
##' data frame.
##'
##' @aliases coef.bayesCox
##' @param object An object returned by function \code{bayesCox}.
##' @param ... Optional arguments. Currently, the only applicable arguemnt is
##'     \code{level} for the credible level. The default value is 0.95.
##' @return A data.frame with 6 columns \code{("Low", "Mid", "High", "Time",
##' "Cov", "Model")}, where \code{"Low"} and \code{"High"} are the posterior
##' 2.5\% and 97.5\% quantiles as default; \code{"Mid"} is the posterior mean;
##' \code{"Cov"} and \code{"Model"} contain character values of the covariates
##' and model types.
##' @seealso \code{\link{bayesCox}}, and \code{\link{plotCoef}}.
##' @keywords extract bayesCox coefficient
##' @examples
##' ## See the examples in bayesCox.
##'
##' @importFrom stats quantile
##' @importFrom data.table .SD
##'
##' @export
coef.bayesCox <- function(object, ...)
{
    ## nonsense to suppress checking notes
    coef <- covariate <- time <- . <- .SD <- NULL

    ## match call for credible level specified
    mcall <- match.call()
    mmcall <- match("level", names(mcall), 0L)
    mcall <- mcall[c(1L, mmcall)]
    mcall[[1L]] <- quote(getLevel)
    level <- eval(mcall)

    ## Monte Carlo samples
    ms <- object$mcmc
    if (is.null(ms)) {
        ms <- read_bayesCox(out = object$out,
                            burn = object$gibbs$burn,
                            thin = object$gibbs$thin)
    }
    betaDat <- bc_beta(ms, object$grid, object$model, object$cov.names)

    if (object$model == "TimeIndep") {
        foo <- function(x) {
            list(time = c(0, object$grid),
                 Low = quantile(x, probs = 0.5 - level / 2, names = FALSE),
                 Mid = mean(x),
                 High = quantile(x, probs = 0.5 + level / 2, names = FALSE),
                 Model = object$model)
        }
        out <- betaDat[, foo(coef), by = covariate]
    } else {
        foo <- function(x) {
            list(Low = quantile(x, probs = 0.5 - level / 2, names = FALSE),
                 Mid = mean(x),
                 High = quantile(x, probs = 0.5 + level / 2, names = FALSE),
                 Model = object$model)
        }
        first_row <- function(x) { x[which.min(time), ] }
        out <- betaDat[, foo(coef), by = list(covariate, time)]
        head_dat <- out[, first_row(.SD), by = covariate]
        head_dat$time <- 0
        out <- rbind(out, head_dat)
        out <- out[order(covariate, time), ]
    }
    data.table::setnames(out, c("covariate", "time"), c("Cov", "Time"))
    out <- out[, c("Low", "Mid", "High", "Time", "Cov", "Model")]
    ## make sure the Cov retains the original order
    out$Cov <- factor(out$Cov, levels = object$cov.names)
    ## return
    attr(out, "level") <- level
    out
}


##' Extract Coefficients from Time-varying Transformation Model
##'
##' Extract coefficient values from \code{tvTran} fitting results, and
##' summarize the point estimate and 95\% credible band into a data frame.
##'
##' @aliases coef.tvTran
##' @param object An object returned by function \code{tvTran}.
##' @param ... Optional arguments. Currently, the only applicable arguemnt is
##'     \code{level} for the credible level. The default value is 0.95.
##' @return A data.frame with 6 columns \code{("Low", "Mid", "High", "Time",
##'     "Cov", "Model")}, where \code{"Mid"} is the point estimates;
##'     \code{"Low"} and \code{"High"} are the 2.5\% and 97.5\% quantiles
##'     estimates from resampling method as default; \code{"Cov"} and
##'     \code{"Model"} contain character values of the covariates and model
##'     type.
##' @seealso \code{\link{tvTran}}, and \code{\link{plotCoef}}.
##' @keywords extract tvTran coefficient
##' @examples
##' ## See the examples in tvTran.
##' @export
coef.tvTran <- function(object, ...)
{
    ## match call for credible level specified
    mcall <- match.call()
    mmcall <- match("level", names(mcall), 0L)
    mcall <- mcall[c(1L, mmcall)]
    mcall[[1L]] <- quote(getLevel)
    level <- eval(mcall)

    K <- object$K
    nBeta <- object$nBeta

    rsMat <- object$rsEst[, seq(1, nBeta * K)]
    betaMat <- cbind(apply(rsMat, 2, quantile, probs = 0.5 - level / 2,
                           na.rm = TRUE, names = FALSE),
                     object$pEst[seq(1, nBeta * K)],
                     apply(rsMat, 2, quantile, probs = 0.5 + level / 2,
                           na.rm = TRUE, names = FALSE))
    ## betaMat[betaMat < -bound | betaMat > bound] <- NA

    ## Insert one more value at time zero
    betaMat <- betaMat[rep(seq(1, nBeta * K),
                           rep(c(2, rep(1, K - 1)), nBeta)), ]

    res <- data.frame(betaMat, rep(c(0, object$eTime), nBeta),
                      rep(object$cov.names, each = K + 1),
                      rep("tvTran", nBeta * (K + 1)))
    colnames(res) <- c("Low", "Mid", "High", "Time", "Cov", "Model")

    ## Make sure the Cov retains the original orde
    res$Cov <- factor(res$Cov, levels = as.character(unique(res$Cov)))
    attr(res, "level") <- level
    res
}


##' Extract Coefficients from Spline Base Cox Model
##'
##' Extract coefficient values from \code{splineCox} fitting results, and
##' summarize the point estimate and 95\% confidence band into a data frame.
##'
##' @aliases coef.splineCox
##' @param object An object returned by function \code{splineCox}.
##' @param ... Optional arguments. Currently, the only applicable arguemnt is
##'     \code{level} for the credible level. The default value is 0.95.
##' @return A data.frame with 6 columns \code{("Low", "Mid", "High", "Time",
##'     "Cov", "Model")}, where \code{"Mid"} is the point estimates;
##'     \code{"Low"} and \code{"High"} are the point estimates plus and minus
##'     1.96 times standard deviations (under default level); \code{"Cov"} and
##'     \code{"Model"} contain character values of the covariates and model
##'     type.
##' @note It essentially expand the break points, and then call function
##'     \code{coxph} in package \code{survival}
##' @seealso \code{\link{splineCox}}, and \code{\link{plotCoef}}.
##' @keywords extract splineCox coefficient
##'
##' @examples
##' ## See the examples in splineCox.
##'
##' @importFrom stats qnorm
##'
##' @export
coef.splineCox <- function(object, ...)
{
    ## match call for credible level specified
    mcall <- match.call()
    mmcall <- match("level", names(mcall), 0L)
    mcall <- mcall[c(1L, mmcall)]
    mcall[[1L]] <- quote(getLevel)
    level <- eval(mcall)

    fit <- object$coxph.fit
    basis <- object$bsp.basis
    K <- 101

    x <- seq(basis$Boundary.knots[1L], basis$Boundary.knots[2L], length = K)
    bspMat <- do.call(splines2::bSpline, c(list(x = x), basis))

    curInd <- 1
    res <- data.frame()
    for (j in seq_along(object$nBeta)) {
        if (!object$is.tv[j]) {
            yMid <- rep(fit$coef[curInd], K)
            ySE <- sqrt(fit$var[curInd, curInd])
            curInd <- curInd + 1
        } else {
            sq <- seq.int(curInd, curInd + basis$df - 1)
            yMid <- c(bspMat %*% fit$coef[sq])
            yVar <- diag(bspMat %*% fit$var[sq, sq] %*% t(bspMat))
            yVar[which(yVar < 0)] <- 0
            ySE <- sqrt(yVar)
            curInd <- curInd + basis$df
        }

        criValue <- qnorm(0.5 + level / 2)
        yLow <- yMid - criValue * ySE
        yHigh <- yMid + criValue * ySE

        res <- rbind(
            res,
            data.frame(Low = yLow, Mid = yMid, High = yHigh,
                       Time = x, Cov = object$cov.names[j],
                       Model = "Spline")
        )
    }

    ## Make sure the Cov retains the original orde
    res$Cov <- factor(res$Cov, levels = as.character(unique(res$Cov)))
    attr(res, "level") <- level
    res
}


### internal function ==========================================================
## help get the possible level specified from ... argument
getLevel <- function(level) {
    if (missing(level))
        return(0.95)
    level
}

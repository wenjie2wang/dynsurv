##
## R package dynsurv by Wenjie Wang, Ming-Hui Chen, Xiaojing Wang, and Jun Yan
## Copyright (C) 2011-2024
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


### Time-varying coefficient transformation model, Peng and Huang (2007) =======
##' Fit Time-varying Transformation Model for Right Censored Survival Data
##'
##' Unlike the time-varying coefficient Cox model, the transformation model
##' fomulates the temporal covariate effects in terms of survival function,
##' i.e., \deqn{S(t|X) = g(\beta_0(t)' X),} where \eqn{g(z) = exp(- exp(z))}.
##' It can be viewed as a functional generalized linear model with response
##' \eqn{I(T > t)}, and other transformation function is possible. The
##' time-varying coefficients are solved a set of estimating equations
##' sequentially.
##'
##' Note that because the time-varying coefficient function is connected to the
##' survival function, it has a different interpretation of the time-varying
##' coefficient function in Cox model.
##'
##' The \code{control} argument is a list of components:
##' \describe{
##'     \item{resample}{A logical value, default \code{TRUE}. If
##'         \code{TRUE}, the model will estimate a 95\% confidence band by
##'         resampling method.}
##'     \item{R}{Number of resamplings, default 30.}
##' }
##' @usage tvTran(formula, data, control = list())
##' @param formula A formula object, with the response on the left of a '~'
##' operator, and the terms on the right. The response must be a survival
##' object as returned by the \code{Surv} function.
##' @param data A data.frame in which to interpret the variables named in the
##' \code{formula}.
##' @param control List of control options.
##' @return An object of S3 class \code{tvTran} representing the fit.
##' @seealso \code{\link{coef.tvTran}}, \code{\link{plotCoef}}.
##'
##' @references
##'
##' Peng, L. and Huang, Y. (2007). Survival analysis with temporal covariate
##' effects. \emph{Biometrika} 94(3), 719--733.
##'
##' @keywords transformation right censor
##' @examples
##' \dontrun{
##' ## Attach the veteran data from the survival package
##' mydata <- survival::veteran
##' mydata$celltype <- relevel(mydata$celltype, ref = "large")
##' myformula <- Surv(time, status) ~ karno + celltype
##'
##' ## Fit the time-varying transformation model
##' fit <- tvTran(myformula, mydata, control = list(resample = TRUE, R = 30))
##'
##' ## Plot the time-varying coefficient function between two time points
##' plotCoef(subset(coef(fit), Time > 15 & Time < 175))
##' }
##' @importFrom stats model.frame model.matrix rnorm
##' @export
tvTran <- function(formula, data, control = list()) {

    Call <- match.call()
    control <- do.call("control_tfun", control)

    ## Right censored data
    mf <- model.frame(formula, data)
    rsp <- mf[, 1]

    ## Event subject index
    eIndex <- which(rsp[, "status"] == 1)

    ## Unique event time
    eTime <- sort(unique(rsp[eIndex, "time"]))
    K <- length(eTime)

    ## Design matrix with intercept term
    X <- model.matrix(formula, data)
    N <- nrow(X)
    nBeta <- ncol(X)
    cov.names <- c("intercept", colnames(X)[-1])
    X <- matrix(X, N, nBeta)

    ## Prepare event, at-risk and offset matrix
    dNMat <- matrix(0, N, K)
    dNMat[eIndex, ] <- outer(rsp[eIndex, "time"], eTime, "==") + 0

    YMat <- outer(rsp[, "time"], eTime, ">=") + 0

    offSetMat <- matrix(0, K, nBeta)

    ## Point estimate
    pEst <- tvTran_lite(X, dNMat, YMat, offSetMat)

    ## Resampling
    if (control$resample) {
        R <- control$R
        rsEst <- matrix(0, R, (nBeta + 1) * K)

        ## Indicator matrix: I(time_i <=  eTime_j) * I(status_i == 1)
        indMat <- outer(rsp[, "time"], eTime, "<=") * rsp[, "status"]

        for (r in 1:R) {
            zeta <- rnorm(N)
            offSetMat <- apply(cbind(0, t(X) %*% (indMat * zeta)), 1, diff)
            rsEst[r, ] <- tvTran_lite(X, dNMat, YMat, offSetMat)
        }
    }
    else
        reEst <- NULL

    ## Return list
    rl <- list(call = Call, eTime = eTime, control = control,
               N = N, K = K, nBeta = nBeta, cov.names = cov.names,
               pEst = pEst, rsEst = rsEst)
    class(rl) <- "tvTran"

    rl
}



### Utility functions ==========================================================
control_tfun <- function(resample = TRUE, R = 30) {
    list(resample = resample, R = R)
}

##' @importFrom nleqslv nleqslv
tvTran_lite <- function(X, dNMat, YMat, offSetMat) {

    ## Estimating function
    f <- function(beta, preExpXb, dN, Y, offSet) {
        c(t(X) %*% (dN - Y * (exp(c(X %*% beta)) - preExpXb))) + offSet
    }

    N <- nrow(X)
    nBeta <- ncol(X)
    K <- ncol(dNMat)

    ## Initialization
    b0 <- rep(0, nBeta)
    betaMat <- matrix(0, K, nBeta)
    preExpXb <- rep(0, N)
    termCode <- rep(0, K)

    ## Solve estimating equations sequentially
    for (k in 1:K) {
        dNVec <- dNMat[, k]
        YVec <- YMat[, k]

        ## Solve time-varying coefficients at time t_i
        res <- nleqslv(b0, f, preExpXb = preExpXb, dN = dNMat[, k],
                       Y = YMat[, k], offSet = offSetMat[k, ], xscalm = "auto")
        beta <- res$x
        termCode[k] <- res$termcd
        betaMat[k, ] <- beta
        preExpXb <- exp(c(X %*% beta))
    }

    betaMat[termCode != 1, ] <- NA

    ## Append termination codes to the end of coefficient estimates
    c(c(betaMat), termCode)
}

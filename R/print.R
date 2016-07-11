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



## ##' Print an object.
## ##'
## ##' An S3 class generic function that displays certain object.
## ##'
## ##' For \code{bayesCox} (\code{tvTran}, or \code{splineCox}) object,
## ##' it prints out the brief summary of the fitted model.
## ##'
## ##' @param x An object used to dispatch a method.
## ##' @param ... Other arguments.



### Print bayesCox object ======================================================
## ## ' @rdname print
## ## ' @aliases print.bayesCox
##' @export
print.bayesCox <- function(x, ...) {

    cat("\nCall:\n")
    dput(x$call)

    cat("\nModel:", x$model, "\n")

    cat("\nPrior for baseline lambda:\n")
    print(data.frame(x$base.prior), row.names = FALSE)

    cat("\nPrior for coefficient beta:\n")
    print(data.frame(x$coef.prior), row.names = FALSE)

    cat("\nOptions for Gibbs sampler:\n")
    print(data.frame(x$gibbs), row.names = FALSE)

    cat("\nOptions for general control:\n")
    print(data.frame(x$control), row.names = FALSE)

    if (is.matrix(x$est$beta))
        beta <- x$est$beta
    else
        beta <- outer(rep(1, x$K), x$est$beta)

    cat("\nBayesian point estimates:\n")
    est <- data.frame(paste("(", c(0, head(x$grid, -1)), ", ",
                           x$grid, "]", sep = ""),
                     cbind(log(x$est$lambda), beta, x$est$jump))
    estNames <- c("interval", "logLambda", paste("beta", x$cov.names, sep = "_"))
    if (!is.null(x$est$jump))
        estNames <- c(estNames, paste("jump", x$cov.names, sep = "_"))

    colnames(est) <- estNames
    print(est, digits = max(options()$digits - 4, 3))

    cat("\nBayesian measures of model fitting:\n")
    print(data.frame(x$measure), row.names = FALSE)
}


### Print tvTran object ========================================================
## ##' @rdname print
## ##' @aliases print.tvTran
##' @export
print.tvTran <- function(x, ...) {

    cat("\nCall:\n")
    dput(x$call)

    cat("\nOptions for control:\n")
    print(data.frame(x$control), row.names = FALSE)

    cat("\nCoefficient estimates:\n")
    est <- data.frame(paste("(", c(0, head(x$eTime, -1)), ", ",
                           x$eTime, "]", sep = ""),
                      matrix(x$pEst[seq(1, x$nBeta * x$K)], nrow = x$K))

    colnames(est) <- c("interval", x$cov.names)
    print(est, digits = max(options()$digits - 4, 3))
}


### Print splineCox object =====================================================
## ##' @rdname print
## ##' @aliases print.splineCox
##' @export
print.splineCox <- function(x, ...) {

    cat("\nCall:\n")
    dput(x$call)

    cat("\nB-spline basis parameters:\n")
    print(data.frame(t(unlist(x$bsp.basis))), row.names = FALSE)

    cat("\nFit results for the expanded data returned by coxph:\n")
    print(x$coxph.fit)
    cat("\n")
}

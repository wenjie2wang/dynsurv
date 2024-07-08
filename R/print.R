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


##' @export
print.bayesCox <- function(x, ...)
{
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
    est <- data.frame(paste0("(", c(0, x$grid[- length(x$grid)]), ", ",
                             x$grid, "]"),
                      cbind(log(x$est$lambda), beta, x$est$jump))
    estNames <- c("interval", "logLambda",
                  paste("beta", x$cov.names, sep = "_"))
    if (!is.null(x$est$jump))
        estNames <- c(estNames, paste("jump", x$cov.names, sep = "_"))

    colnames(est) <- estNames
    print(est, digits = max(options()$digits - 4, 3))

    cat("\nBayesian measures of model fitting:\n")
    print(data.frame(x$measure), row.names = FALSE)
    ## return x invisibly
    invisible(x)
}


##' @export
print.tvTran <- function(x, ...) {

    cat("\nCall:\n")
    dput(x$call)

    cat("\nOptions for control:\n")
    print(data.frame(x$control), row.names = FALSE)

    cat("\nCoefficient estimates:\n")
    est <- data.frame(paste0("(", c(0, x$eTime[- length(x$eTime)]), ", ",
                             x$eTime, "]"),
                      matrix(x$pEst[seq_len(x$nBeta * x$K)], nrow = x$K))

    colnames(est) <- c("interval", x$cov.names)
    print(est, digits = max(options()$digits - 4, 3))
    ## return x invisibly
    invisible(x)
}


##' @export
print.splineCox <- function(x, ...) {

    cat("\nCall:\n")
    dput(x$call)

    cat("\nB-spline basis parameters:\n")
    print(data.frame(t(unlist(x$bsp.basis))), row.names = FALSE)

    cat("\nFit results for the expanded data returned by coxph:\n")
    print(x$coxph.fit)
    cat("\n")
    ## return x invisibly
    invisible(x)
}

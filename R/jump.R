##
## R package dynsurv by Wenjie Wang, Ming-Hui Chen, Xiaojing Wang, and Jun Yan
## Copyright (C) 2011-2022
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


##' Generic function for jump information
##'
##' @param object An object returned by function \code{bayesCox}.
##' @param ... Other arguments.
##' @export
jump <- function(object, ...) UseMethod("jump", object)


##' @describeIn jump Extract Jump Information from Bayesian Dyanmic Model
##'
##' Extract number of coefficient pieces from \code{bayesCox} fitting results,
##' and summarize them into a data frame. It is only applicable when
##' \code{model="Dynamic"} is specified.
##'
##' @aliases jump jump.bayesCox
##' @return A data.frame with 3 columns \code{("Count", "Iter", "Cov")}, where
##' \code{"Count"} is the number of coefficient pieces (jumps) for each
##' iteration; \code{Iter} is the iteration number; \code{Cov} contains the
##' character values of the covariates.
##' @seealso \code{\link{bayesCox}}, \code{\link{plotJumpTrace}}, and
##' \code{\link{plotJumpHist}}.
##' @keywords extract bayesCox jump
##' @examples
##' ## See the examples in bayesCox
##'
##' @importFrom stats model.frame
##'
##' @export
jump.bayesCox <- function(object, ...)
{
    ## nonsense to supprese checking notes
    mcmc.sample <- covariate <- NULL

    ## Monte Carlo samples
    ms <- object$mcmc
    if (is.null(ms)) {
        ms <- read_bayesCox(out = object$out,
                            burn = object$gibbs$burn,
                            thin = object$gibbs$thin)
    }
    out <- bc_jump(ms, object$grid, object$model, object$cov.names)
    out <- out[, list(Count = sum(jump)), by = list(mcmc.sample, covariate)]
    data.table::setnames(out, c("mcmc.sample", "covariate"), c("Iter", "Cov"))
    out$Model <- object$model
    out$Cov <- factor(out$Cov, levels = object$cov.names)
    ## return
    out[, c("Iter", "Cov", "Model", "Count")]
}


### Extract the latent variance nu from "bayesCox" object ======================
##' Generic function for the latent variance of coefficients
##'
##' @param object An object returned by function \code{bayesCox}.
##' @param ... Other arguments.
##' @export
nu <- function(object, ...) UseMethod("nu", object)


##' @describeIn nu Extract Latent Variance from Bayesian Cox Model
##'
##' Extract latent variance of coefficients from \code{bayesCox} fitting
##' results, and summarize them into a data frame. It is applicable when
##' \code{model="TimeVarying"} or \code{model="Dynamic"}, and
##' \code{coef.prior=list(type="HAR1")}.
##'
##' For details, see section on prior model in Wang (2013) and Wang (2014).
##' The latent variance of coefficients in prior model was denoted as omega
##' in Wang (2013).
##'
##' @aliases nu.bayesCox
##' @return A data.frame with 4 columns \code{("Iter", "Model", "Cov",
##' "Value")}, where \code{Iter} is the iteration number; \code{Model} and
##' \code{Cov} contain the character values of the model type and covariates.
##' @seealso \code{\link{bayesCox}}, and \code{\link{plotNu}}.
##' @references
##' X. Wang, M.-H. Chen, and J. Yan (2013). Bayesian dynamic regression
##' models for interval censored survival data with application to children
##' dental health. Lifetime data analysis, 19(3), 297--316.
##'
##' X. Wang, X. Sinha, J. Yan, and M.-H. Chen (2014). Bayesian inference of
##' interval-censored survival data. In: D. Chen, J. Sun, and K. Peace,
##' Interval-censored time-to-event data: Methods and applications, 167--195.
##' @examples
##' ## See the examples in bayesCox.
##'
##' @importFrom stats model.frame
##'
##' @export
nu.bayesCox <- function(object, ...)
{
    ## Monte Carlo samples
    ms <- object$mcmc
    if (is.null(ms)) {
        ms <- read_bayesCox(out = object$out,
                            burn = object$gibbs$burn,
                            thin = object$gibbs$thin)
    }
    out <- bc_nu(ms, object$grid, object$model, object$cov.names)
    data.table::setnames(out,
                         c("mcmc.sample", "covariate", "nu"),
                         c("Iter", "Cov", "Value"))
    out$Model <- object$model
    ## return
    out[, c("Iter", "Model", "Cov", "Value")]
}

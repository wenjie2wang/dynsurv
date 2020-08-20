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


##' Get the MCMC Samples from \code{bayesCox}
##'
##' Returns the MCMC samples produced by
##' \code{bayesCox} into data frames.
##'
##' @param object A \code{bayesCox} object
##' @param parts A character vector specifying the parts to be exacted from the
##'     MCMC output text file produced by \code{bayesCox}.  One or more
##'     following options can be specified: \code{"h0"} for baseline hazard
##'     function, \code{"coef"} for covariate coefficients, \code{"nu"} for
##'     sampled latent variance of coefficients, \code{"jump"} for indicators of
##'     jumps, and \code{"all"} for all of the above.  The default value is
##'     \code{c("h0", "beta")}.
##' @param ... Other arguments that are not used now.
##'
##' @importFrom stats reshape
##'
##' @export
bayesCoxMcmc <- function(object,
                         parts = c("h0", "coef"),
                         ...)
{
    if (! is.bayesCox(object)) {
        stop("The input must be a 'bayesCox' object.")
    }
    parts <- match.arg(parts, choices = c("h0", "coef", "nu", "jump", "all"),
                       several.ok = TRUE)
    ## Monte Carlo samples
    ms <- read.table(file = object$out)
    ms <- ms[seq.int(object$gibbs$burn + 1, nrow(ms), by = object$gibbs$thin), ]
    row.names(ms) <- NULL
    iter <- nrow(ms)
    sample_id <- seq_len(iter)
    ## dimension of baseline
    grid_ <- object$grid
    K <- length(grid_)
    cK <- if (object$model == "TimeIndep") 1 else K
    nBeta <- length(object$cov.names)
    ## baseline hazard
    h0Dat <- NULL
    if (any(c("h0", "all") %in% parts)) {
        h0Dat <- ms[, seq_len(K), drop = FALSE]
        h0Dat$mcmc.sample <- sample_id
        h0Dat <- stats::reshape(h0Dat,
                                idvar = "mcmc.sample",
                                varying = seq_len(K),
                                v.names = "h0",
                                times = grid_,
                                direction = "long")
        h0Dat <- h0Dat[order(h0Dat$mcmc.sample), ]
        row.names(h0Dat) <- NULL
        attr(h0Dat, "reshapeLong") <- NULL
        h0Dat <- h0Dat
    }
    ## covariate coefficients
    betaDat <- NULL
    if (any(c("coef", "all") %in% parts)) {
        betaDat <- do.call(rbind, lapply(seq_len(nBeta), function(j) {
            jj <- K + (j - 1) * cK + 1
            if (object$model != "TimeIndep") {
                out <- ms[, seq.int(jj, jj + cK - 1), drop = FALSE]
                out$mcmc.sample <- sample_id
                out <- stats::reshape(out,
                                      idvar = "mcmc.sample",
                                      varying = seq_len(cK),
                                      v.names = "coef",
                                      times = grid_,
                                      direction = "long")
                out <- out[order(out$mcmc.sample), ]
                attr(out, "reshapeLong") <- NULL
            } else {
                out <- ms[, jj, drop = FALSE]
                out$mcmc.sample <- sample_id
                colnames(out)[1L] <- "coef"
                out <- out[, c("mcmc.sample", "coef")]
            }
            out$covariate <- object$cov.names[j]
            row.names(out) <- NULL
            out
        }))
    }
    ## latent variance of coefficients
    nuDat <- NULL
    if (any(c("nu", "all") %in% parts)) {
        ## only for time-varying and dynamic model and coef.prior = HAR1
        if (! object$model %in% c("TimeVarying", "Dynamic")) {
            warning(
                "Latent variance of coefficients is only available for",
                "the time-varying coefficient models or the dynamic models."
            )
        } else {
            nuDat <- do.call(rbind, lapply(seq_len(nBeta), function(j) {
                jj <- (1 + nBeta) * K + j
                out <- ms[, jj, drop = FALSE]
                out$mcmc.sample <- sample_id
                colnames(out)[1L] <- "nu"
                out <- out[, c("mcmc.sample", "nu")]
                out$covariate <- object$cov.names[j]
                row.names(out) <- NULL
                out
            }))
        }
    }
    ## jump
    jumpDat <- NULL
    if (any(c("jump", "all") %in% parts)) {
        ## only for dynamic model
        if (object$model != "Dynamic") {
            warning(
                "Jump indicators are only available for the dynamic models."
            )
        } else {
            jumpDat <- do.call(rbind, lapply(seq_len(nBeta), function(j) {
                jj <- (j + nBeta) * K + nBeta + 1
                out <- ms[, seq.int(jj, jj + K - 1), drop = FALSE]
                out$mcmc.sample <- sample_id
                out <- stats::reshape(out,
                                      idvar = "mcmc.sample",
                                      varying = seq_len(K),
                                      v.names = "jump",
                                      times = grid_,
                                      direction = "long")
                out <- out[order(out$mcmc.sample), ]
                attr(out, "reshapeLong") <- NULL
                out$covariate <- object$cov.names[j]
                row.names(out) <- NULL
                out
            }))
        }
    }
    ## return
    structure(
        list(
            h0 = h0Dat,
            coef = betaDat,
            nu = nuDat,
            jump = jumpDat
        ),
        class = c("bayesCoxMcmc")
    )
}

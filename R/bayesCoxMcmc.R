##
## R package dynsurv by Wenjie Wang, Ming-Hui Chen, Xiaojing Wang, and Jun Yan
## Copyright (C) 2011-2023
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
    ms <- object$mcmc
    if (is.null(ms)) {
        ms <- read_bayesCox(object$out, object$gibbs$burn, object$gibbs$thin)
    }
    ## baseline hazard
    h0Dat <- NULL
    if (any(c("h0", "all") %in% parts)) {
        h0Dat <- bc_h0(ms, object$grid, object$model, object$cov.names)
    }
    ## covariate coefficients
    betaDat <- NULL
    if (any(c("coef", "all") %in% parts)) {
        betaDat <- bc_beta(ms, object$grid, object$model, object$cov.names)
    }
    ## latent variance of coefficients
    nuDat <- NULL
    if (any(c("nu", "all") %in% parts)) {
        ## only for time-varying and dynamic model and coef.prior = HAR1
        if (! object$model %in% c("TimeVarying", "Dynamic")) {
            warning(
                "Latent variance of coefficients is only available for ",
                "the time-varying coefficient models or the dynamic models."
            )
        } else {
            nuDat <- bc_nu(ms, object$grid, object$model, object$cov.names)
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
            jumpDat <- bc_jump(ms, object$grid, object$model, object$cov.names)
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


### internal functions =======================================================

## read mcmc outputs
##' @importFrom data.table fread
read_bayesCox <- function(out, burn, thin)
{
    ms <- fread(
        file = out,
        sep = " ",
        showProgress = FALSE,
        verbose = FALSE
    )
    ms <- ms[seq.int(burn + 1, nrow(ms), by = thin), ]
    row.names(ms) <- NULL
    ms
}

## tidy mcmc
##' @importFrom data.table melt data.table setnames rbindlist

## baseline hazard function
bc_h0 <- function(ms, grid, model, cov.names)
{
    ## set dimension
    iter <- nrow(ms)
    sample_id <- seq_len(iter)
    K <- length(grid)
    cK <- if (model == "TimeIndep") 1 else K
    nBeta <- length(cov.names)
    ## wide to long
    h0Dat <- ms[, seq_len(K), drop = FALSE, with = FALSE]
    h0Dat$mcmc.sample <- sample_id
    h0Dat <- data.table::melt(h0Dat,
                              idvars = "mcmc.sample",
                              measure.vars = seq_len(K),
                              variable.name = "tmp",
                              value.name = "h0",
                              variable.factor = FALSE)
    tmp <- data.table(
        tmp = paste0("V", seq_len(K)),
        time = grid
    )
    h0Dat <- merge(h0Dat, tmp, by = "tmp")
    h0Dat$tmp <- NULL
    ## return
    h0Dat[, c("mcmc.sample", "time", "h0")]
}

## covariate coefficients
bc_beta <- function(ms, grid, model, cov.names)
{
    ## set dimension
    iter <- nrow(ms)
    sample_id <- seq_len(iter)
    K <- length(grid)
    cK <- if (model == "TimeIndep") 1 else K
    nBeta <- length(cov.names)
    ## wide to long
    if (model == "TimeIndep") {
        idx <- seq.int(K + 1, K + nBeta)
        betaDat <- ms[, idx, drop = FALSE, with = FALSE]
        betaDat$mcmc.sample <- sample_id
        betaDat <- data.table::melt(betaDat,
                                    idvars = "mcmc.sample",
                                    measure.vars = seq_len(nBeta),
                                    variable.name = "covariate",
                                    value.name = "coef",
                                    variable.factor = TRUE)
        betaDat$covariate <- factor(
            betaDat$covariate,
            levels = levels(betaDat$covariate),
            labels = cov.names
        )
        return(betaDat)
    }
    betaDat <- rbindlist(lapply(seq_len(nBeta), function(j) {
        jj <- K + (j - 1) * cK + 1
        idx <- seq.int(jj, jj + cK - 1)
        out <- ms[, idx, drop = FALSE, with = FALSE]
        out$mcmc.sample <- sample_id
        out <- data.table::melt(out,
                                idvars = "mcmc.sample",
                                measure.vars = seq_len(K),
                                variable.name = "tmp",
                                value.name = "coef",
                                variable.factor = TRUE)
        tmp <- data.table(
            tmp = paste0("V", idx),
            time = grid
        )
        out <- merge(out, tmp, by = "tmp")
        out$tmp <- NULL
        out$covariate <- cov.names[j]
        row.names(out) <- NULL
        out
    }))
    betaDat$covariate <- factor(betaDat$covariate, levels = cov.names)
    ## return
    betaDat[, c("mcmc.sample", "time", "coef", "covariate")]
}

## latent variance
bc_nu <- function(ms, grid, model, cov.names)
{
    nuDat <- NULL
    if (model %in% c("TimeVarying", "Dynamic")) {
        ## set dimension
        iter <- nrow(ms)
        sample_id <- seq_len(iter)
        K <- length(grid)
        nBeta <- length(cov.names)
        ## wide to long
        nuDat <- ms[, seq.int((1 + nBeta) * K + 1, (1 + nBeta) * K + nBeta),
                    with = FALSE]
        nuDat$mcmc.sample <- sample_id
        nuDat <- melt(nuDat,
                      idvars = "mcmc.sample",
                      measure.vars = seq_len(nBeta),
                      variable.name = "covariate",
                      value.name = "nu",
                      variable.factor = TRUE)
        nuDat$covariate <- factor(
            nuDat$covariate,
            levels = levels(nuDat$covariate),
            labels = cov.names
        )
    }
    ## return
    nuDat
}

## jump
bc_jump <- function(ms, grid, model, cov.names)
{
    jumpDat <- NULL
    if (model == "Dynamic") {
        ## set dimension
        iter <- nrow(ms)
        sample_id <- seq_len(iter)
        K <- length(grid)
        nBeta <- length(cov.names)
        ## wide to long
        jumpDat <- rbindlist(lapply(seq_len(nBeta), function(j) {
            jj <- (j + nBeta) * K + nBeta + 1
            idx <- seq.int(jj, jj + K - 1)
            out <- ms[, idx, drop = FALSE, with = FALSE]
            out$mcmc.sample <- sample_id
            tmp <- data.table(
                tmp = paste0("V", idx),
                time = grid
            )
            out <- melt(out,
                        idvars = "mcmc.sample",
                        measure.vars = seq_len(K),
                        variable.name = "tmp",
                        value.name = "jump",
                        variable.factor = TRUE)

            out <- merge(out, tmp, by = "tmp")
            out$tmp <- NULL
            out$covariate <- cov.names[j]
            row.names(out) <- NULL
            out
        }))
    }
    jumpDat$covariate <- factor(jumpDat$covariate, levels = cov.names)
    ## return
    jumpDat[, c("mcmc.sample", "time", "jump", "covariate")]
}

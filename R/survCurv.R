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


##' Estimated Survival Function or Cumulative Hazard Function
##'
##' Estimated survival function or cumulative hazard function from posterior
##' sample for an object returned by function \code{bayesCox}.
##'
##' The estimated survival curve is a step function representing the posterior
##' mean survival proportion at the given time grid from the posterior sample.
##' The credible interval for the survival curve is constructed based on the
##' quantiles of all the survival curves from posterior sample at given credible
##' level. More details were available in Section posterior computation of Wang
##' (2016).
##'
##' @aliases survCurve
##'
##' @param object An object returned by function \code{bayesCox}.
##' @param newdata An optional data frame used to generate a design matrix.
##' @param type An optional character value indicating the type of function to
##'     compute. The possible values are "survival" and "cumhaz". The former
##'     means the estimated survival function; the latter represents the
##'     estimated cumulative hazard function for the given \code{newdata}.
##' @param level A numerical value between 0 and 1 indicating the level of
##'     cradible band.
##' @param centered A logical value. If \code{TRUE}, the mean function for the
##'     given \code{newdata} will be computed. The default is \code{FALSE}.
##' @param ... Other arguments for further usage.
##'
##' @return A data frame with column: "Low", "Mid", "High", "Time", "Design",
##'     and "type", and attribute, "surv" valued as "survCurve".
##'
##' @seealso
##' \code{\link{bayesCox}},
##' \code{\link{survDiff}}, and
##' \code{\link{plotSurv}}.
##'
##' @references
##'
##' Wang, W., Chen, M. H., Chiou, S. H., Lai, H. C., Wang, X., Yan, J.,
##' & Zhang, Z. (2016). Onset of persistent pseudomonas aeruginosa infection in
##' children with cystic fibrosis with interval censored data.
##' \emph{BMC Medical Research Methodology}, 16(1), 122.
##'
##' @examples
##' ## See the examples in bayesCox.
##'
##' @importFrom stats model.frame model.matrix delete.response terms setNames
##'
##' @export
survCurve <- function(object, newdata, type = c("survival", "cumhaz"),
                      level = 0.95, centered = FALSE, ...)
{
    ## nonsense, just to suppress Note from R CMD check --as-cran
    `(Intercept)` <- NULL

    ## generate design matrix from newdata
    tt <- stats::terms(object$formula)
    Terms <- stats::delete.response(tt)
    nBeta <- object$nBeta
    if (missing(newdata)) {
        X <- matrix(0, nrow = 1, ncol = nBeta)
        colnames(X) <- object[["cov.names"]]
        if (object$control$intercept)
            X[, "intercept"] <- 1
        rownames(X) <- "Baseline"
    } else {
        mf <- stats::model.frame(Terms, newdata, xlev = object$xlevels)
        X <- stats::model.matrix(Terms, mf)
        ## remove intercept and deplicated rows
        X <- unique(base::subset(X, select = - `(Intercept)`))
        if (object$control$intercept)
            X <- cbind(intercept = 1, X)
        if (ncol(X) != nBeta) {
            stop("The number of input covariates does not ",
                 "match with the 'bayesCox' object")
        }
    }
    nDesign <- nrow(X)

    ## read ms
    ms <- object$mcmc
    if (is.null(ms)) {
        ms <- read_bayesCox(out = object$out,
                            burn = object$gibbs$burn,
                            thin = object$gibbs$thin)
    }
    ms <- as.matrix(ms)

    iter <- nrow(ms)
    bGrid <- c(0, object$grid)
    K <- length(object$grid)
    deltaT <- diff(bGrid)
    h0Mat <- ms[, seq_len(K), drop = FALSE]
    cK <- ifelse(object$model == "TimeIndep", 1, K)
    betaMat <- as.matrix(ms[, seq(K + 1, K + nBeta * cK)])

    ## function to generate ht and Ht for each design
    compHt <- function(oneX) {
        oneX <- matrix(oneX, nrow = nBeta)
        expXbeta <- matrix(0, ncol = cK, nrow = iter)
        ## for j_th time point on the grid
        for (j in seq_len(cK)) {
            betat <- betaMat[, seq(j, by = cK, length.out = nBeta)]
            expXbeta[, j] <- exp(rowMeans(betat %*% oneX))
        }
        ht <- h0Mat * as.numeric(expXbeta)
        Ht <- matrix(0, ncol = K + 1L, nrow = iter)
        for (i in seq_len(iter)) {
            Ht[i, - 1L] <- cumsum(ht[i, ] * deltaT)
        }
        Ht
    }

    if (centered) {
        Ht <- compHt(oneX = t(X))
        nDesign <- 1
    } else {
        Ht <- matrix(NA, nrow = iter, ncol = (K + 1) * nDesign)
        for (i in seq(nDesign)) {
            Ht[, seq(1 + (i - 1) * (K + 1), i * (K + 1))] <-
                compHt(oneX = X[i, ])
        }
    }

    ## function type
    type <- match.arg(type)
    if (type == "cumhaz") {
        ## cumulative hazard function
        HtQT <- data.frame(t(apply(Ht, 2, ciBand, level = level)))
        colnames(HtQT) <- c("Low", "Mid", "High")
        HtQT$Time <- rep(bGrid, times = nDesign)
        HtQT$Design <- if (centered)
                           "centered"
                       else
                           factor(rep(rownames(X), each = K + 1))

        HtQT$type <- "cumhaz"
        attr(HtQT, "surv") <- "survCurve"
        return(HtQT)
    }
    ## else survival function
    St <- exp(- Ht)
    StQT <- data.frame(t(apply(St, 2, ciBand, level = level)))
    colnames(StQT) <- c("Low", "Mid", "High")
    StQT$Time <- rep(bGrid, times = nDesign)
    StQT$Design <- if (centered)
                       "centered"
                   else
                       factor(rep(rownames(X), each = K + 1))

    StQT$type <- "survival"
    attr(StQT, "surv") <- "survCurve"
    ## return
    StQT
}


##' Estimated Difference Between Survival or Cumulative Hazard Functions
##'
##' \code{survDiff} returns estimated survival function or cumulative function
##' from posterior estimates. Note that currently, the function is only
##' applicable to the Bayesian dynamic Cox model with dynamic hazard, where the
##' control argument is specified to be \code{control = list(intercept = TRUE)}
##' in function \code{bayesCox}.
##'
##' The estimated difference between survival curves is a step function
##' representing the difference between the posterior mean survival proportion
##' at the given time grid from the posterior sample.  Its credible interval is
##' constructed based on the quantiles of all the pair difference between the
##' survival curves from posterior sample at given credible level.
##'
##' @aliases survDiff
##'
##' @param object An object returned by function \code{bayesCox}.
##' @param newdata An optional data frame used to generate a design matrix.
##'     Note that it must lead to a design matrix with two different design.
##' @param type An optional character value indicating the type of function to
##'     compute. The possible values are "survival" and "cumhaz". The former
##'     means the estimated survival function; the latter represents the
##'     estimated cumulative hazard function for the given \code{newdata}.
##' @param level A numerical value between 0 and 1 indicating the level of
##'     cradible band.
##'
##' @param ... Other arguments for further usage.
##' @return A data frame with column: "Low", "Mid", "High", "Time", "Design",
##'     and "type", and attribute, "surv" valued as "survDiff".
##'
##' @references
##'
##' Wang, W., Chen, M. H., Chiou, S. H., Lai, H. C., Wang, X., Yan, J.,
##' & Zhang, Z. (2016). Onset of persistent pseudomonas aeruginosa infection in
##' children with cystic fibrosis with interval censored data.
##' \emph{BMC Medical Research Methodology}, 16(1), 122.
##'
##' @seealso
##' \code{\link{bayesCox}}, \code{\link{survCurve}}, and \code{\link{plotSurv}}.
##'
##' @examples
##' ## See the examples in bayesCox.
##'
##' @importFrom stats model.frame model.matrix delete.response terms
##'
##' @export
survDiff <- function(object, newdata, type = c("survival", "cumhaz"),
                     level = 0.95, ...)
{
    ## nonsense, just to suppress Note from R CMD check --as-cran
    `(Intercept)` <- NULL

    ## generate design matrix from newdata
    tt <- stats::terms(object$formula)
    Terms <- stats::delete.response(tt)
    nBeta <- object$nBeta
    if (missing(newdata)) {
        stop("The argument 'newdata' is missing.")
    } else {
        mf <- stats::model.frame(Terms, newdata, xlev = object$xlevels)
        X <- stats::model.matrix(Terms, mf)
        ## remove intercept and deplicated rows
        X <- unique(base::subset(X, select = -`(Intercept)`))
        if (object$control$intercept)
            X <- cbind(intercept = 1, X)
        if (ncol(X) != nBeta) {
            stop("The number of input covariates does not ",
                 "match with the 'bayesCox' object")
        }
    }
    nDesign <- nrow(X)
    if (nDesign != 2) {
        stop("The 'newdata' must contain two different designs.")
    }

    ## read ms
    ms <- object$mcmc
    if (is.null(ms)) {
        ms <- read_bayesCox(out = object$out,
                            burn = object$gibbs$burn,
                            thin = object$gibbs$thin)
    }
    ms <- as.matrix(ms)

    iter <- nrow(ms)
    bGrid <- c(0, object$grid)
    K <- length(object$grid)
    deltaT <- diff(bGrid)
    h0Mat <- ms[, seq_len(K), drop = FALSE]
    cK <- ifelse(object$model == "TimeIndep", 1, K)
    nBeta <- length(object$cov.names)
    betaMat <- as.matrix(ms[, seq(K + 1, K + nBeta * cK)])

    ## function to generate ht and Ht for each design
    compHt <- function(oneX) {
        oneX <- matrix(oneX, nrow = nBeta)
        expXbeta <- matrix(0, ncol = cK, nrow = iter)
        ## for j_th time point on the grid
        for (j in seq_len(cK)) {
            betat <- betaMat[, seq(j, by = cK, length.out = nBeta)]
            expXbeta[, j] <- exp(rowMeans(betat %*% oneX))
        }
        ht <- h0Mat * as.numeric(expXbeta)
        Ht <- matrix(0, ncol = K + 1L, nrow = iter)
        for (i in seq_len(iter)) {
            Ht[i, - 1L] <- cumsum(ht[i, ] * deltaT)
        }
        Ht
    }

    res1 <- compHt(oneX = X[1L, ])
    res2 <- compHt(oneX = X[2L, ])
    Designdiff <- paste(rownames(X)[2L], rownames(X)[1L], sep = " vs. ")

    ## function type
    type <- match.arg(type)

    ## cumulative hazard function difference
    if (type == "cumhaz") {
        Htdiff <- res2 - res1
        HtdiffQT <- data.frame(t(apply(Htdiff, 2, ciBand, level = level)))
        colnames(HtdiffQT) <- c("Low", "Mid", "High")
        HtdiffQT$Time <- bGrid
        HtdiffQT$Design <- Designdiff
        HtdiffQT$type <- "cumhaz"
        attr(HtdiffQT, "surv") <- "survDiff"
        return(HtdiffQT)
    }

    ## survival function difference
    St1 <- exp(- res1)
    St2 <- exp(- res2)
    Stdiff <- St2 - St1
    StdiffQT <- data.frame(t(apply(Stdiff, 2, ciBand, level = level)))
    colnames(StdiffQT) <- c("Low", "Mid", "High")
    StdiffQT$Time <- bGrid
    StdiffQT$Design <- Designdiff
    StdiffQT$type <- "survival"
    attr(StdiffQT, "surv") <- "survDiff"
    ## else return
    StdiffQT
}


##' Plot Survival Curves (or Cumulative Hazard Function) and their difference
##'
##' Plot the survival curves (or cumulative hazard) and their difference for
##' objects returned by function \code{survCurve} or \code{survDiff}.  By using
##' \code{ggplot2} plotting system, the plots generated are able to be further
##' customized properly.
##'
##' @aliases plotSurv
##'
##' @param object An object returned by function \code{survCurve} or
##'     \code{survDiff}.
##' @param legendName An optional name for the figure legend.
##' @param conf.int A logical value indicating whether to plot the credible
##'     interval(s).
##' @param smooth A logical value, default \code{FALSE}. If \code{TRUE}, plot
##'     the coefficients as smooth lines; otherwise, plot the coefficients as
##'     piece-wise constant step functions.
##' @param lty An optional numeric vector indicating line types specified to
##'     different groups: 0 = blank, 1 = solid, 2 = dashed, 3 = dotted, 4 =
##'     dotdash, 5 = longdash, 6 = twodash.
##' @param col An optional character or numeric vector indicating line colors
##'     specified to different groups.
##' @param ... Other arguments for future usage.
##'
##' @return A \code{ggplot} object.
##'
##' @seealso
##' \code{\link{bayesCox}}, \code{\link{survCurve}}, and
##' \code{\link{survDiff}}.
##'
##' @examples
##' ## See the examples in bayesCox.
##'
##' @importFrom ggplot2 ggplot aes_string aes geom_step scale_color_manual
##'     scale_linetype_manual geom_line ylab ggtitle
##'
##' @export
plotSurv <- function(object, legendName = "", conf.int = FALSE,
                     smooth = FALSE, lty, col, ...)
{
    ## nonsense, just to suppress Note from R CMD check --as-cran
    Mid <- Low <- High <- Time <- NULL

    ## about linetypes
    ## 0 = blank, 1 = solid, 2 = dashed, 3 = dotted,
    ## 4 = dotdash, 5 = longdash, 6 = twodash.

    ## eliminate possible empty levels
    Design <- factor(object$Design)
    nDesign <- length(levels(Design))

    ## set line types and colors
    lty <- if (missing(lty))
               setNames(rep(1, nDesign), levels(Design))
           else
               setNames(lty[seq(nDesign)], levels(Design))

    col <- if (missing(col))
               gg_color_hue(nDesign)
           else
               col[seq(nDesign)]

    ## initialize ggplot object
    p <- ggplot(data = object, aes_string(x = "Time"))

    if (! smooth) {
        if (nDesign == 1) {
            p <-  p + geom_step(mapping = aes(x = Time, y = Mid))

            if (conf.int)
                p <- p +
                    geom_step(mapping = aes(x = Time, y = Low),
                              linetype = "3313") +
                    geom_step(mapping = aes(x = Time, y = High),
                              linetype = "3313")
        } else {
            p <- p + geom_step(mapping = aes(x = Time, y = Mid, color = Design,
                                             linetype = Design)) +
                scale_color_manual(values = col, name = legendName) +
                scale_linetype_manual(values = lty, name = legendName)

            if (conf.int)
                p <- p +
                    geom_step(mapping = aes(x = Time, y = Low, color = Design),
                              linetype = "3313") +
                    geom_step(mapping = aes(x = Time, y = High, color = Design),
                              linetype = "3313")
        }
    } else {
        if (nDesign == 1) {
            p <-  p + geom_line(mapping = aes(x = Time, y = Mid))

            if (conf.int)
                p <- p +
                    geom_line(mapping = aes(x = Time, y = Low),
                              linetype = "3313") +
                    geom_line(mapping = aes(x = Time, y = High),
                              linetype = "3313")
        } else {
            p <- p + geom_line(mapping = aes(x = Time, y = Mid, color = Design,
                                             linetype = Design)) +
                scale_color_manual(values = col, name = legendName) +
                scale_linetype_manual(values = lty, name = legendName)

            if (conf.int)
                p <- p +
                    geom_line(mapping = aes(x = Time, y = Low, color = Design),
                              linetype = "3313") +
                    geom_line(mapping = aes(x = Time, y = High, color = Design),
                              linetype = "3313")
        }
    }
    type <- unique(object$type)[1L]
    if (attr(object, "surv") == "survDiff") {
        p <- if (type == "survival")
                 p + ylab("Estimated Survival Proportion Difference") +
                     ggtitle(unique(Design))
             else
                 p + ylab("Estimated Cumulative Hazard Difference") +
                     ggtitle(unique(Design))
    } else if (attr(object, "surv") == "survCurve") {
        p <- if (type == "survival")
                 p + ylab("Estimated Survival Proportion")
             else
                 p + ylab("Estimated Cumulative Hazard")
    } else
        stop("Unknown 'surv' attribute.")
    p
}


### internal functions =========================================================
## function to compute ci and posterior mean
##' @importFrom stats quantile
ciBand <- function(x, level = 0.95) {
    lev <- (1 - level) / 2
    c(quantile(x, probs = lev, names = FALSE), mean(x),
      quantile(x, probs = 1 - lev, names = FALSE))
}


## function to emulate the default colors used in ggplot2
##' @importFrom grDevices hcl
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[seq_len(n)]
}

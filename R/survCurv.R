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



### Generate the surviva curves for dynamic model with dynamic hazard.
### read out file into R (possibly the most time-consuming part)

##' Estimated Survival Function or Cumulative Hazard Function
##'
##' \code{survCurve} returns estimated suvival function or cumulative function
##' from posterior sample. Note that the function is currently only
##' applicable to the Bayesian dynamic Cox model with dynamic hazard, where the
##' control argument is specified to be \code{control = list(intercept = TRUE)}
##' in function \code{bayesCox}.
##'
##' The estimated survival curve is a step function representing the posterior
##' mean survival proportion at the given time grid from the posterior sample.
##' The credible interval for the survival curve is constructed based on the
##' quantiles of all the survival curves from posterior sample at given credible
##' level. More details were available in Section posterior computation of
##' Wang (2016).
##'
##' @aliases survCurve
##' @usage
##' survCurve(object, newdata, type = c("survival", "cumhaz"),
##'           level = 0.95, centered = FALSE, cache = FALSE, ...)
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
##' @param cache A logical value. If \code{TRUE}, the cache RData file will be
##'     generated in the working directory to improve the performance of
##'     function \code{survCurve} and \code{survDiff}. This option would be
##'     quite helpful if the number of MCMC sample is large.
##' @param ... Other arguments for further usage.
##' @return A data frame with column: "Low", "Mid", "High", "Time", "Design",
##' and "type", and attribute, "surv" valued as "survCurve".
##' @seealso
##' \code{\link{bayesCox}},
##' \code{\link{survDiff}}, and
##' \code{\link{plotSurv}}.
##' @references
##' Wang, W., Chen, M. H., Chiou, S. H., Lai, H. C., Wang, X., Yan, J.,
##' & Zhang, Z. (2016). Onset of persistent pseudomonas aeruginosa infection in
##' children with cystic fibrosis with interval censored data.
##' \emph{BMC Medical Research Methodology}, 16(1), 122.
##' @examples
##' ## See the examples in bayesCox.
##' @importFrom stats model.frame model.matrix delete.response terms setNames
##' @importFrom utils read.table
##' @export
survCurve <- function(object, newdata, type = c("survival", "cumhaz"),
                      level = 0.95, centered = FALSE, cache = FALSE, ...){

    ## NOTE:
    ##
    ## The construction of crediable band fails for time-varying model
    ## and dynamic model with time-varying model at reference level.
    ## This is bacause covariates at reference level are all set to be 0.
    ## And the posterior sample of lambda is not included in outputs
    ## from function bayesCox.
    ##
    ## It may be fixed later. Then this function can be expanded to other
    ## two models completely.

    ## ONLY applicable to Bayesian dynamic Cox model with dynamic hazard
    if (! (object$control$intercept & object$model == "Dynamic")) {
        stop(paste("'survCurve' is only applicable to",
                   "Bayesian dynamic Cox model with dynamic hazard."))
    }

    ## nonsense, just to suppress Note from R CMD check --as-cran
    `(Intercept)` <- NULL

    ## generate design matrix from newdata
    tt <- stats::terms(object$formula)
    Terms <- stats::delete.response(tt)
    nBeta <- object$nBeta
    if (missing(newdata)) {
        X <- matrix(0, nrow = 1, ncol = nBeta)
        colnames(X) <- object[["cov.names"]]
        ## if (object$control$intercept)
        X[, "intercept"] <- 1
        rownames(X) <- "Baseline"
    } else {
        mf <- stats::model.frame(Terms, newdata, xlev = object$xlevels)
        X <- stats::model.matrix(Terms, mf)
        ## remove intercept and deplicated rows
        X <- unique(base::subset(X, select = -`(Intercept)`))
        ## if (object$control$intercept)
        X <- cbind(intercept = 1, X)
        if (ncol(X) != nBeta) {
            stop(paste("The number of input covariates does not",
                       "match with the 'bayesCox' object"))
        }
    }
    nDesign <- nrow(X)

    ## check or generate cache and ms
    cacheName <- paste(substring(object$out, 1, nchar(object$out) - 4),
                       "cache", sep = "_")
    cacheOut <- paste(cacheName, ".RData", sep = "")
    if (! file.exists(cacheOut)) {
        ms <- as.matrix(read.table(file = object$out))
        dimnames(ms) <- NULL
        ms <- ms[seq(object$gibbs$burn + 1, nrow(ms), by = object$gibbs$thin), ]
        assign(cacheName, ms)
        if (cache) save(list = cacheName, file = cacheOut)
    } else {
        load(cacheOut)
        if (! cacheName %in% ls()) {
            stop("Cache file does not match. Please remove local cache files.")
        }
        assign("ms", get(cacheName))
    }

    iter <- nrow(ms)
    bGrid <- c(0, object$grid)
    K <- length(bGrid) - 1
    deltaT <- diff(bGrid)
    cK <- ifelse(object$model == "TimeIndep", 1, K)
    betaMat <- as.matrix(ms[, seq(K + 1, K + nBeta * cK)])

    ## function to generate ht and Ht for each design
    foo <- function(oneX){
        oneX <- matrix(oneX, nrow = nBeta)
        tempHt <- tempht <- matrix(0, ncol = cK + 1, nrow = iter)
        ## for j_th time point on the grid
        for (j in seq(2, cK + 1)) {
            betat <- betaMat[, seq(j - 1, by = cK, length.out = nBeta)]
            tempht[, j] <- object$est$lambda[j - 1] *
                exp(rowMeans(betat %*% oneX))
            tempHt[, j] <- as.matrix(tempht[, 2 : j]) %*% deltaT[seq(j - 1)]
        }
        tempHt
    }
    if (centered) {
        Ht <- foo(oneX = t(X))
        nDesign <- 1
    } else {
        Ht <- matrix(NA, nrow = iter, ncol = (cK + 1) * nDesign)
        for (i in seq(nDesign)) {
            Ht[, seq(1 + (i - 1) * (cK + 1), i * (cK + 1))] <-
                foo(oneX = X[i, ])
        }
    }

    ## function type
    type <- match.arg(type)
    if (type == "cumhaz"){
        ## cumulative hazard function
        HtQT <- data.frame(t(apply(Ht, 2, ciBand, level = level)))
        colnames(HtQT) <- c("Low", "Mid", "High")
        HtQT$Time <- rep(bGrid, times = nDesign)
        HtQT$Design <- if (centered) {
                           "centered"
                       } else {
                           factor(rep(rownames(X), each = cK + 1))
                       }
        HtQT$type <- "cumhaz"
        attr(HtQT, "surv") <- "survCurve"
        return(HtQT)
    }
    ## else survival function
    St <- exp(- Ht)
    StQT <- data.frame(t(apply(St, 2, ciBand, level = level)))
    colnames(StQT) <- c("Low", "Mid", "High")
    StQT$Time <- rep(bGrid, times = nDesign)
    StQT$Design <- if (centered) {
                       "centered"
                   } else {
                       factor(rep(rownames(X), each = cK + 1))
                   }
    StQT$type <- "survival"
    attr(StQT, "surv") <- "survCurve"
    ## return
    StQT
}



##' Estimated Difference Between Survival or Cumulative Hazard Functions
##'
##' \code{survDiff} returns estimated suvival function or cumulative function
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
##' @usage survDiff(object, newdata, type = c("survival", "cumhaz"),
##'          level = 0.95, cache = FALSE, ...)
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
##' @param cache A logical value. If \code{TRUE}, the cache RData file will be
##'     generated in the working directory to improve the performance of
##'     function \code{survCurve} and \code{survDiff}. This option would be
##'     quite helpful if the number of MCMC sample is large.
##' @param ... Other arguments for further usage.
##' @return A data frame with column: "Low", "Mid", "High", "Time", "Design",
##'     and "type", and attribute, "surv" valued as "survDiff".
##' @seealso \code{\link{bayesCox}}, \code{\link{survCurve}}, and
##' \code{\link{plotSurv}}.
##' @examples
##' ## See the examples in bayesCox.
##' @importFrom stats model.frame model.matrix delete.response terms
##' @importFrom utils read.table
##' @export
survDiff <- function (object, newdata, type = c("survival", "cumhaz"),
                      level = 0.95, cache = FALSE, ...) {

    ## ONLY applicable to Bayesian dynamic Cox model with dynamic hazard
    if (! (object$control$intercept & object$model == "Dynamic")) {
        stop(paste("'survCurve' is only applicable to",
                   "Bayesian dynamic Cox model with dynamic hazard."))
    }

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
        ## if (object$control$intercept)
        X <- cbind(intercept = 1, X)
        if (ncol(X) != nBeta) {
            stop(paste("The number of input covariates does not",
                       "match with the 'bayesCox' object"))
        }
    }
    nDesign <- nrow(X)
    if (nDesign != 2) {
        stop("The 'newdata' must contain two different designs.")
    }

    ## check or generate cache and ms
    cacheName <- paste(substring(object$out, 1, nchar(object$out) - 4),
                       "cache", sep = "_")
    cacheOut <- paste(cacheName, ".RData", sep = "")
    if (! file.exists(cacheOut)) {
        ms <- as.matrix(read.table(file = object$out))
        dimnames(ms) <- NULL
        ms <- ms[seq(object$gibbs$burn + 1, nrow(ms), by = object$gibbs$thin), ]
        assign(cacheName, ms)
        if (cache) save(list = cacheName, file = cacheOut)
    } else {
        load(cacheOut)
        if (! cacheName %in% ls()) {
            stop("Cache file does not match. Please remove local cache files.")
        }
        assign("ms", get(cacheName))
    }

    iter <- nrow(ms)
    bGrid <- c(0, object$grid)
    K <- length(bGrid) - 1
    deltaT <- diff(bGrid)
    cK <- ifelse(object$model == "TimeIndep", 1, K)
    nBeta <- length(object$cov.names)
    betaMat <- as.matrix(ms[, seq(K + 1, K + nBeta * cK)])
    ## function to generate ht and Ht for each design
    foo <- function(oneX){
        oneX <- matrix(oneX, nrow = nBeta)
        tempHt <- tempht <- matrix(0, ncol = cK + 1, nrow = iter)
        ## for j_th time point on the grid
        for(j in seq(2, cK + 1)){
            betat <- betaMat[, seq(j - 1, by = cK, length.out = nBeta)]
            tempht[, j] <- object$est$lambda[j - 1] *
                exp(rowMeans(betat %*% oneX))
            tempHt[, j] <- as.matrix(tempht[, 2 : j]) %*% deltaT[seq(j - 1)]
        }
        tempHt
    }
    res1 <- foo(oneX = X[1, ])
    res2 <- foo(oneX = X[2, ])
    Designdiff <- paste(rownames(X)[2], rownames(X)[1], sep=" vs. ")

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
##' @usage plotSurv(object, legendName = "", conf.int = FALSE, smooth = FALSE,
##'          lty, col, ...)
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
##' @return A \code{ggplot} object.
##' @seealso \code{\link{bayesCox}}, \code{\link{survCurve}}, and
##' \code{\link{survDiff}}.
##' @examples
##' ## See the examples in bayesCox.
##' @importFrom ggplot2 ggplot aes_string aes geom_step scale_color_manual
##'     scale_linetype_manual geom_line ylab ggtitle
##' @export
plotSurv <- function(object, legendName = "", conf.int = FALSE,
                     smooth = FALSE, lty, col, ...){

    ## nonsense, just to suppress Note from R CMD check --as-cran
    Mid <- Low <- High <- Time <- NULL

    ## about linetypes
    ## 0 = blank, 1 = solid, 2 = dashed, 3 = dotted,
    ## 4 = dotdash, 5 = longdash, 6 = twodash.

    ## eliminate possible empty levels
    Design <- factor(object$Design)
    nDesign <- length(levels(Design))

    ## set line types and colors
    lty <- if (missing(lty)) {
               setNames(rep(1, nDesign), levels(Design))
           } else {
               setNames(lty[seq(nDesign)], levels(Design))
           }
    col <- if (missing(col)) {
               gg_color_hue(nDesign)
           } else {
               col[seq(nDesign)]
           }

    ## initialize ggplot object
    p <- ggplot(data = object, aes_string(x = "Time"))

    if (! smooth) {
        if (nDesign == 1) {
            p <-  p + geom_step(mapping = aes(x = Time, y = Mid))
            if (conf.int) {
                p <- p +
                    geom_step(mapping = aes(x = Time, y = Low),
                              linetype = "3313") +
                    geom_step(mapping = aes(x = Time, y = High),
                              linetype = "3313")
            }
        } else {
            p <- p + geom_step(mapping = aes(x = Time, y = Mid, color = Design,
                                             linetype = Design)) +
                scale_color_manual(values = col, name = legendName) +
                scale_linetype_manual(values = lty, name = legendName)
            if (conf.int) {
                p <- p +
                    geom_step(mapping = aes(x = Time, y = Low, color = Design),
                              linetype = "3313") +
                    geom_step(mapping = aes(x = Time, y = High, color = Design),
                              linetype = "3313")
            }
        }
    } else {
        if (nDesign == 1) {
            p <-  p + geom_line(mapping = aes(x = Time, y = Mid))
            if (conf.int) {
                p <- p +
                    geom_line(mapping = aes(x = Time, y = Low),
                              linetype = "3313") +
                    geom_line(mapping = aes(x = Time, y = High),
                              linetype = "3313")
            }
        } else {
            p <- p + geom_line(mapping = aes(x = Time, y = Mid, color = Design,
                                             linetype = Design)) +
                scale_color_manual(values = col, name = legendName) +
                scale_linetype_manual(values = lty, name = legendName)
            if (conf.int) {
                p <- p +
                    geom_line(mapping = aes(x = Time, y = Low, color = Design),
                              linetype = "3313") +
                    geom_line(mapping = aes(x = Time, y = High, color = Design),
                              linetype = "3313")
            }
        }
    }
    if (attr(object, "surv") == "survDiff") {
        if (all(unique(object$type) == "survival")) {
            p <- p + ylab("Estimated Survival Proportion Difference") +
                ggtitle(unique(Design))
        } else if (all(unique(object$type) == "cumhaz")) {
            p <- p + ylab("Estimated Cumulative Hazard Difference") +
                ggtitle(unique(Design))
        }
    } else if (attr(object, "surv") == "survCurve") {
        if (all(unique(object$type) == "survival")) {
            p <- p + ylab("Estimated Survival Proportion")
        } else if (all(unique(object$type) == "cumhaz")) {
            p <- p + ylab("Estimated Cumulative Hazard")
        }
    }
    p
}



### internal functions =========================================================
## function to compute ci and posterior mean
##' @importFrom stats quantile
ciBand <- function(x, level = 0.95){
    c(quantile(x, probs = (1 - level) / 2, names = FALSE), mean(x),
      quantile(x, probs = 1 - (1 - level) / 2, names = FALSE))
}


## function to emulate the default colors used in ggplot2
##' @importFrom grDevices hcl
gg_color_hue <- function(n){
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1 : n]
}

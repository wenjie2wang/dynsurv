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


### Fit a time-varying coefficient Cox model, using B-splines ==================
##' Fit Spline Based Cox Model for Right Censored Survival Data
##'
##' Rearrange the rignt censored survival data in a counting process style.
##' Model the time-varying coefficient function using B-splines. The fit is
##' done by introducing pseudo time-dependent covariates and then calling
##' function \code{coxph} in \pkg{survival} package.
##'
##' The \code{control} argument is a list of components:
##' \describe{
##'     \item{df:}{
##'         degree of freedom for the B-splines, default 5;}
##'     \item{knots:}{interior knots point, default \code{NULL}. If
##'         \code{NULL}, the knots will be automatically choosen;}
##'     \item{boundary:}{lower and upper boundaries for the spline
##'         function, default \code{NULL}. If \code{NULL}, the minimun
##'         and maximun finite event time or censoring time will be
##'         specified.}
##' }
##'
##' @param formula A formula object, with the response on the left of a '~'
##' operator, and the terms on the right. The response must be a survival
##' object as returned by the \code{Surv} function.
##' @param data A data.frame in which to interpret the variables named in the
##' \code{formula}.
##' @param control List of control options.
##' @return An object of S3 class \code{splineCox} representing the fit.
##' @note This function is essentially a wrapper function of \code{coxph} for
##' the expanded data set. It does not implements the algorithm disscussed in
##' the reference paper. These authors implemented their algorithm into a
##' \code{tvcox} package, which is more efficient for larger data set, but may
##' not be stable compared to \code{coxph}.
##' @seealso \code{\link{coef.splineCox}}, \code{\link{plotCoef}}.
##'
##' @references
##'
##' Perperoglou, A., le Cessie, S., & van Houwelingen, H. C. (2006). A fast
##' routine for fitting Cox models with time varying effects of the
##' covariates. Computer Methods and Programs in Biomedicine, 81(2), 154--161.
##'
##' @keywords B-spline Cox right censor
##'
##' @examples
##' \dontrun{
##' ## Attach the veteran data from the survival package
##' mydata <- survival::veteran
##' mydata$celltype <- relevel(mydata$celltype, ref = "large")
##' myformula <- Surv(time, status) ~ karno + celltype
##'
##' ## Fit the time-varying transformation model
##' fit <- splineCox(myformula, mydata, control = list(df = 5))
##'
##' ## Plot the time-varying coefficient function between two time points
##' plotCoef(subset(coef(fit), Time > 15 & Time < 175), smooth = TRUE)
##' }
##' @importFrom stats model.matrix model.frame as.formula
##' @importFrom survival coxph
##' @importFrom splines2 bSpline
##' @export
splineCox <- function(formula, data, control = list()) {

    Call <- match.call()
    control <- do.call("control_sfun", control)

    ## Model matrix after expanding the factor covariate
    mm <- model.matrix(formula, data)
    N <- nrow(mm)
    nBeta <- ncol(mm) - 1
    is.tv <- rep(TRUE, nBeta)
    is.tv[grep("const[(]([A-Za-z0-9._]*)[)]", colnames(mm)[-1])] <- FALSE

    cov.names <- gsub("const[(]([A-Za-z0-9._]*)[)]", "\\1", colnames(mm)[-1])

    ## Model frame before expanding the factor covariate
    mf <- model.frame(formula, data)
    nCov <- ncol(mf) - 1

    FtNms <- rep("Ft:", nCov)
    FtNms[grep("const[(]([A-Za-z0-9._]*)[)]", names(mf)[-1])] <- ""

    names(mf) <- gsub("const[(]([A-Za-z0-9._]*)[)]", "\\1", names(mf))

    ## First 3 columns  =  c("id", "time", "status")
    DF <- cbind(id = 1:N, mf[, 1][, 1:2], mf[, -1, drop = FALSE])

    ## Prepare B-spline paramters
    boundary <- control$boundary
    if (is.null(boundary))
        boundary <- range(DF$time)

    knots <- control$knots
    df <- control$df
    if (is.null(control$knots)) {
        insideTime <- subset(DF$time, DF$time >=  boundary[1] &
                                      DF$time <=  boundary[2])

        ## number of interior knots  =  df - degree(3) - intercept(1)
        sq <- seq.int(from = 0, to = 1, length.out = df - 2)[- c(1, df-2)]
        knots <- stats::quantile(insideTime, sq)
    }

    basis <- list(df = control$df, knots = knots, intercept = TRUE,
                  Boundary.knots = boundary)

    ## Call coxph to fit the expanded data
    newDF <- expand(DF, id = "id", time = "time", status = "status")

    ## B-spline basis matrix
    Ft <- do.call(splines2::bSpline, c(list(x = newDF$tStop), basis))

    newFml <- as.formula(paste("survival::Surv(tStart, tStop, status) ~ ",
                               paste(paste(FtNms, names(mf)[-1], sep = ""),
                                     collapse = "+"), "+ cluster(id)"))

    fit <- coxph(newFml, newDF)

    rl <- list(call = Call, control = control, bsp.basis = basis,
               N = N, nBeta = nBeta, cov.names = cov.names, is.tv = is.tv,
               coxph.fit = fit)

    class(rl) <- "splineCox"
    rl
}


### Utility functions ==========================================================
##' @importFrom data.table rbindlist
expand <- function(data, id = "id", time = "time", status = "status")
{
    pos <- match(c(id, time, status), names(data))
    if (length(pos) !=  3)
        stop("Variable names not match!\n")

    eventTime <- sort(unique(data[, time][data[, status] == 1]))

    foo <- function(x) {
        tStop <- union(subset(eventTime, eventTime <= max(x[, time])),
                       max(x[, time]))
        tStart <- c(0, tStop[- length(tStop)])
        st <- rep(0, length(tStop))
        st[tStop %in% x[x[, status] == 1, time]] <- 1

        cbind(x[1, -pos[2:3]], tStart = tStart, tStop = tStop,
              st = st, row.names = NULL)
    }

    res <- data.table::rbindlist(by(data, data$id, foo))
    names(res)[ncol(res)] <- status
    res
}

control_sfun <- function(df = 5, knots = NULL, boundary = NULL)
{
    list(df = df, knots = knots, boundary = boundary)
}

const <- function(x) x

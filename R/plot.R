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

##############################################################################
# Plot coefficient data frame
##############################################################################


##' Plot Coefficient Function
##' 
##' Plot coefficient values formatted in a data frame returned by function
##' \code{coef}.
##' 
##' To plot estimated coefficient functions from different models together, one
##' can first combine the data frames returned by \code{coef}, and then call
##' \code{plotCoef}, for example, \code{plotCoef(rbind(coef(fit1),
##' coef(fit2)))}.
##' 
##' To specify the time range of the plot, one can either utilize the
##' \code{ggplot} functionality, say \code{plotCoef(coef(fit)) + xlim(2, 10)};
##' or manipulate the data frame first before calling \code{plotCoef}, e.g.,
##' \code{plotCoef(subset(coef(fit), Time > 2 & Time < 10))}.
##' 
##' @usage plotCoef(object, smooth=FALSE, \dots{})
##' @param object a data.frame returned by function \code{coef}.
##' @param smooth a logical value, default \code{FALSE}. If \code{TRUE}, plot
##' the coefficients as smooth lines; otherwise, plot the coefficients as
##' piece-wise constant step functions.
##' @param \dots other arguments.
##' @return A \code{ggplot} object.
##' @seealso \code{\link{coef.bayesCox}}, \code{\link{coef.splineCox}}, and
##' \code{\link{coef.tvTran}}.
##' @keywords plot coefficient
##' @examples
##' 
##' # See the examples in bayesCox, splineCox, and tvTran.
##' 
##' @export plotCoef
plotCoef <- function(object, smooth=FALSE, ...) {

    p <- ggplot(data=object, aes_string(x="Time"))

    if (!smooth)
        p <- p + geom_step(aes_string(y="Mid"), direction="vh") +
            geom_step(aes_string(y="High"), direction="vh", linetype=2) +
            geom_step(aes_string(y="Low"), direction="vh", linetype=2)
    else
       p <- p + geom_line(aes_string(y="Mid")) +
           geom_line(aes_string(y="High"), linetype=2) +
           geom_line(aes_string(y="Low"), linetype=2)

    if (length(levels(factor(object$Model))) == 1)
        p <- p + facet_wrap(~ Cov, scales="free_y")
    else
        p <- p + facet_grid(Cov ~ Model, scales="free_y")

    p <- p + ylab("Coefficient") +
        theme(plot.margin=unit(rep(0, 4), "lines"))

    p
}


##############################################################################
# Plot iteration jump data frame returned by jump.bayesCox
##############################################################################
plotJumpTrace <- function(object, ...) {
    p <- ggplot(data=object, aes_string(x="Iter", y="Count")) +
        geom_line(size=0.1, alpha=0.6) +
        facet_wrap(~ Cov) +
        xlab("Iteration") + ylab("Pieces of Coefficient") +
        theme(plot.margin=unit(rep(0, 4), "lines"))

    p
}

plotJumpHist <- function(object, ...) {
  ## p <- ggplot(data=object, aes(x=factor(Count))) +
  p <- ggplot(data=object, aes_string(x="Count")) +
    ##   stat_bin(aes(y =..count../sum(..count..))) +
    stat_bin(aes_string(y = "..density..")) +
      facet_wrap(~ Cov) +
        xlab("Pieces of Coefficient") + ylab("Relative Frequency") +
          theme(plot.margin=unit(rep(0, 4), "lines"))
  
  p
}

##############################################################################
# Plot the latent variance nu from the bayesCox model
##############################################################################


##' Plot Latent Variance in Bayesian Cox Model
##' 
##' Plot the latent variance \code{nu} when the hierarchical AR(1) process
##' prior is used for the \code{bayesCox} model. It is applicable when
##' \code{model="TimeVarying"} or \code{model="Dynamic"}, and
##' \code{coef.prior=list(type="HAR1")}. The input data frame is returned by
##' function \code{nu}.
##' 
##' 
##' @usage plotNu(object, \dots{})
##' @param object a data.frame returned by the function \code{nu}.
##' @param \dots other arguments.
##' @return A \code{ggplot} object.
##' @seealso \code{\link{nu.bayesCox}}.
##' @keywords plot latent variance
##' @examples
##' 
##' # See the examples in bayesCox
##' 
##' @export plotNu
plotNu <- function(object, ...) {
  cnt <- "..density.."
  ## p <- ggplot(data=object, aes_string(x="Value")) +
  p <- ggplot(data=object, aes_string(x="Value")) +
    ##  stat_bin(aes(y =..count../sum(..count..))) +
    stat_bin(aes_string(y = cnt)) +
      xlab("Nu") + ylab("Relative Frequency") +
        theme(plot.margin=unit(rep(0, 4), "lines"))
  
  if (length(levels(factor(object$Model))) == 1)
    p <- p + facet_wrap(~ Cov)
  else
    p <- p + facet_grid(Cov ~ Model)

  p
}

##############################################################################
##
##   R package dynsurv by Xiaojing Wang, Jun Yan, and Ming-Hui Chen
##   Copyright (C) 2011
##
##   This file is part of the R package dynsurv.
##
##   The R package dynsurv is free software: you can redistribute it and/or
##   modify it under the terms of the GNU General Public License as published
##   by the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package dynsurv is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package dynsurv. If not, see <http://www.gnu.org/licenses/>.
##
##############################################################################

##############################################################################
# Extract jump data frame from "bayesCox" object
##############################################################################
# setGeneric("jump", function(object, ...) standardGeneric("jump"))
jump <- function(object, ...) UseMethod("jump", object)



##' Extract Jump Information from Bayesian Dyanmic Model
##' 
##' Extract number of coefficient pieces from \code{bayesCox} fitting results,
##' and summarize them into a data frame. It is only applicable when
##' \code{model="Dynamic"} is specified.
##' 
##' 
##' @aliases jump jump.bayesCox
##' @usage \method{jumpbayesCox}(object, \dots{})
##' @param object an object returned by function \code{bayesCox}.
##' @param \dots other arguments.
##' @return A data.frame with 3 columns \code{("Count", "Iter", "Cov")}, where
##' \code{"Count"} is the number of coefficient pieces (jumps) for each
##' iteration; \code{Iter} is the iteration number; \code{Cov} contains the
##' character values of the covariates.
##' @seealso \code{\link{bayesCox}}, \code{\link{plotJumpTrace}}, and
##' \code{\link{plotJumpHist}}.
##' @keywords extract bayesCox jump
##' @examples
##' 
##' # See the examples in bayesCox
##' 
##' @export jump.bayesCox
jump.bayesCox <- function(object, ...) {
    # Monte Carlo samples
    ms <- as.matrix(read.table(file=object$out))
    dimnames(ms) <- NULL
    ms <- ms[seq(object$gibbs$burn + 1, nrow(ms), by=object$gibbs$thin), ]
    iter <- nrow(ms)

    grid <- object$grid
    K <- length(grid)
    covNms <- object$cov.names
    nBeta <- length(covNms)

    jumpMat <- as.matrix(ms[, seq((1 + nBeta) * K + nBeta + 1, (1 + 2 * nBeta) * K + nBeta)])

    # Number of jumps for each iteration
    csMat <- diag(1, nBeta, nBeta)
    csMat <- as.matrix(csMat[rep(1:nBeta, each=K), ])
    iterJump <- data.frame(c(jumpMat %*% csMat), rep(1:iter, nBeta), rep(covNms, each=iter))
    colnames(iterJump) <- c("Count", "Iter", "Cov")
    iterJump$Cov <- factor(iterJump$Cov, levels=as.character(unique(iterJump$Cov)))

    # Number of jumps at each time grid point
    # timeJump <- data.frame(colMeans(jumpMat), rep(grid, nBeta), rep(covNms, each=K))
    # colnames(timeJump) <- c("Jump", "Time", "Cov")
    # timeJump <- timeJump[-seq(K, nBeta * K, by=K), ]
    # list(iterJump=iterJump, timeJump=timeJump)

    iterJump
}

##############################################################################
# Extract the latent variance nu from "bayesCox" object
##############################################################################
# setGeneric("nu", function(object, ...) standardGeneric("nu"))
nu <- function(object, ...) UseMethod("nu", object)



##' Extract Latent Variance from Bayesian Cox Model
##' 
##' Extract latent variance from \code{bayesCox} fitting results, and summarize
##' them into a data frame. It is applicable when \code{model="TimeVarying"} or
##' \code{model="Dynamic"}, and \code{coef.prior=list(type="HAR1")}.
##' 
##' 
##' @aliases nu nu.bayesCox
##' @usage \method{nubayesCox}(object, \dots{})
##' @param object an object returned by function \code{bayesCox}.
##' @param \dots other arguments.
##' @return A data.frame with 4 columns \code{("Iter", "Model", "Cov",
##' "Value")}, where \code{Iter} is the iteration number; \code{Model} and
##' \code{Cov} contain the character values of the model type and covariates.
##' @seealso \code{\link{bayesCox}}, and \code{\link{plotNu}}.
##' @keywords extract bayesCox latent variance
##' @examples
##' 
##' # See the examples in bayesCox
##' 
##' @export nu.bayesCox
nu.bayesCox <- function(object, ...) {
    # Monte Carlo samples
    ms <- as.matrix(read.table(file=object$out))
    dimnames(ms) <- NULL
    ms <- ms[seq(object$gibbs$burn + 1, nrow(ms), by=object$gibbs$thin), ]
    iter <- nrow(ms)

    K <- length(object$grid)
    covNms <- object$cov.names
    nBeta <- length(covNms)

    nuMat <- as.matrix(ms[, seq((1 + nBeta) * K + 1, (1 + nBeta) * K + nBeta)])
    res <- data.frame(1:iter, object$model, nuMat)
    colnames(res) <- c("Iter", "Model", covNms)

    res <- melt(res, c("Iter", "Model"))
    colnames(res) <- c("Iter", "Model", "Cov", "Value")
    res
}

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


##' Breast Cancer Data
##'
##' The breast cancer data from Finkelstein (1985) has been analyzed
##' extensively for illustrating new methods in modeling interval censored
##' data. The objective of the study was to compare the time to cosmetic
##' deterioration between two groups: 46 patients receiving radiotherapy only
##' and 48 patients receiving radiotherapy plus chemotherapy. Because the
##' detection of deterioration required a clinic visit, the 56 women who
##' experience deterioration were interval censored, and the other38 women who
##' did not were right censored.
##' @name bcos
##' @aliases bcos bcos.grid
##' @docType data
##' @format \code{bcos} is a data frame with 94 observations and 3 columns.
##' \describe{
##'     \item{left:}{left censoring time;}
##'     \item{right:}{right censoring time;}
##'     \item{trt:}{treatment (\code{Rad} = radiotherapy only,
##'          \code{RadChem} = radiotherapy plus chemotherapy).}
##' }
##'
##' \code{bcos.grid} is a numeric vector of grid time points.
##' @references D.M. Finkelstein, and R.A. Wolfe (1985). A semiparametric model
##' for regression analysis of interval-censored failure time data.
##' \emph{Biometrics} 41: 731--740.
##'
##' @keywords datasets
##' @examples
##' data(bcos)
NULL


##' Tooth Data
##'
##' The tooth data was from a longitudinal prospective dental study performed
##' in Flanders (Belgium) in 1996 -- 2001. Every one of 4,386 randomly sampled
##' children in the cohort was examined annually by one of 16 trained dentists,
##' resulting at most 6 dental observations for each child. The outcome of
##' interest was the time to emergence of permanent tooth 24, which was either
##' interval censored (2,775, 63\%) or right censored (1,611, 37\%).
##'
##' @name tooth
##' @aliases tooth tooth.grid
##' @docType data
##' @format \code{tooth} is a data frame with 4,386 observations and 7 columns
##' \describe{
##'     \item{id:}{children's id;}
##'     \item{left:}{left censoring time;}
##'     \item{right:}{right censoring time where infinity is coded as 999.}
##'     \item{sex:}{gender of children (0 = boy, 1 = girl);}
##'     \item{dmf:}{status of the primary predecessor of this tooth (0 = sound,
##'         1 = delayed, missing or filled);}
##'     \item{rightInf:}{right censoring time where infinity is coded as
##'         \code{Inf};}
##'     \item{rightNA:}{right censoring time where infinity is coded as
##'         \code{NA}.}
##' }
##' \code{tooth.grid} is a numeric vector of grid time points.
##' @references J. Vanobbergen, L. Martens, D. Declerck, and M. Lesaffre
##' (2000). The signal tandmobiel(r) project: a longitudinal intervention oral
##' health promotion study in Flanders (Belgium): baseline and first year
##' results. \emph{European Journal of Paediatric Dentistry} 2, 87.
##'
##' G. Gomez, M. Calle, R. Oller, and K. Langohr (2009). Tutorial on methods
##' for interval-censored data and their implementation in R. \emph{Statistical
##' Modeling} 9(4), 259.
##' @source Adapted from the data set available at
##' \url{https://grass.upc.edu/en/software/tooth24}.
##' @keywords datasets
##' @examples
##' data(tooth)
NULL

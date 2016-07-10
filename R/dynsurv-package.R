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


#' dynsurv: Time-varying coefficient models for interval censored data.
#'
#' Functions to fit time-varying coefficient models for interval censored and
#' right censored survival data. Three major approaches are implemented:
#' \enumerate{
#' \item Bayesian Cox model with time-independent, time-varying or dynamic
#'     coefficients for right censored and interval censored data;
#' \item Spline based time-varying coefficient Cox model for right censored
#'     data;
#' \item Transformation model with time-varying coefficients for right censored
#'     data using estimating equations.
#' }
#'
#' @docType package
#' @name dynsurv
#' @importFrom survival Surv
#' @useDynLib dynsurv
NULL

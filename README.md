# R Package dynsurv

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/dynsurv)](http://cran.r-project.org/package=dynsurv)

The R package **dynsurv** provides functions to fit time-varying coefficient
models for interval censored and right censored survival data.

Three major approaches are implemented:

* 1) Bayesian Cox model with
time-independent, time-varying or dynamic coefficients for
right censored and interval censored data;

* 2) Spline based time-varying coefficient Cox model for right censored data;

* 3) Transformation model with time-varying coefficients for
right censored data using estimating equations.


## Installation

You may install the released version from
[CRAN](http://cran.rstudio.com/package=dynsurv):

```r
install.packages("dynsurv", dependencies = TRUE)
```


## Development

[![Build Status](https://travis-ci.org/wenjie2wang/dynsurv.svg?branch=dev)](https://travis-ci.org/wenjie2wang/dynsurv)

The latest version of package is under development
at [GitHub](https://github.com/wenjie2wang/dynsurv) in branch 'dev'.  You may
consider installing the latest version with the help of **devtools** if it is
able to pass the building check by Travis CI.

```R
if (! require(devtools)) install.packages("devtools", dependencies = TRUE)
devtools::install_git("git://github.com/wenjie2wang/dynsurv.git", branch = "dev")
```


## Usage

See [package help manual](https://cran.rstudio.com/web/packages/dynsurv/dynsurv.pdf)
for details and demonstration.


## License

The R package reda is free software: You can redistribute it and/or
modify it under the terms of the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License, or
any later version (at your option).
See the [GNU General Public License](http://www.gnu.org/licenses/) for details.

The R package reda is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

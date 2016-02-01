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

You may install the stable version on
[CRAN](http://cran.rstudio.com/package=dynsurv):

```r
install.packages("dynsurv")
```


## Usage

```r
help(pacakge = "dynsurv")
library(dynsurv)
## help on main function for model fitting
?bayesCox
```

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

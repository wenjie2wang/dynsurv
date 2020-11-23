# dynsurv

[![CRAN_Status_Badge][r-pkg-badge]][cran-url]
[![Build Status][gha-icon]][gha-url]

The R package **dynsurv** provides functions fitting time-varying coefficient
models for interval censored and right censored survival data.

Three major approaches are implemented:

1. Bayesian Cox model with time-independent, time-varying or dynamic
   coefficients for right censored and interval censored data;
1. Spline based time-varying coefficient Cox model for right censored data;
1. Transformation model with time-varying coefficients for right censored data
   using estimating equations.


## Installation

You may install the released version from [CRAN][cran-url].

```R
install.packages("dynsurv")
```


## Development

The latest version of package is under development at [GitHub][github-url].  If
it is able to pass the automated package checks, one may install it by

```R
if (! require(remotes)) install.packages("remotes")
remotes::install_github("wenjie2wang/dynsurv", upgrade = "never")
```

## Get Started

```
library(dynsurv)
?bayesCox
```

## License

[GNU General Public License][gpl] (≥ 3)


[r-pkg-badge]: https://www.r-pkg.org/badges/version/dynsurv
[cran-url]: https://CRAN.R-project.org/package=dynsurv
[github-url]: https://github.com/wenjie2wang/dynsurv
[gha-icon]: https://github.com/wenjie2wang/dynsurv/workflows/R-CMD-check/badge.svg
[gha-url]: https://github.com/wenjie2wang/dynsurv/actions
[codecov]: https://codecov.io/gh/wenjie2wang/dynsurv
[codecov-main]: https://codecov.io/gh/wenjie2wang/dynsurv/branch/main/graph/badge.svg
[gpl]: https://www.gnu.org/licenses/

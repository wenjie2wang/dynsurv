\dontrun{

library(dynsurv)
set.seed(1216)

### load one of the following two data sets

## breast cancer data
data(bcos)       # load bcos and bcos.grid
mydata <- bcos
mygrid <- bcos.grid
myformula <- Surv(left, right, type = "interval2") ~ trt

## tooth data
## data(tooth)   # load tooth and tooth.grid
## mydata <- tooth
## mygrid <- tooth.grid
## myformula <- Surv(left, rightInf, type = "interval2") ~ dmf + sex

### Fit time-independent coefficient model
fit0 <- bayesCox(myformula, mydata, out = "tiCox.txt", model = "TimeIndep",
                 base.prior = list(type = "Gamma", shape = 0.1, rate = 0.1),
                 coef.prior = list(type = "Normal", mean = 0, sd = 1),
                 gibbs = list(iter = 100, burn = 20, thin = 1,
                              verbose = TRUE, nReport = 20))
## plot coefficient estimates
plotCoef(coef(fit0, level = 0.9))

## Plot the estimated survival function for given new data
newDat <- data.frame(trt = c("Rad", "RadChem"))
row.names(newDat) <- unique(newDat$trt)
plotSurv(survCurve(fit0, newDat), conf.int = TRUE)

## optional: get MCMC samples for more estimates of interest
mcmc_list <- read.bayesCox(fit0)
str(mcmc_list)

### Fit time-varying coefficient model
fit1 <- bayesCox(myformula, mydata, out = "tvCox.txt", model = "TimeVary",
                 base.prior = list(type = "Gamma", shape = 0.1, rate = 0.1),
                 coef.prior = list(type = "AR1", sd = 1),
                 gibbs = list(iter = 100, burn = 20, thin = 1,
                              verbose = TRUE, nReport = 20))
plotCoef(coef(fit1))
plotSurv(survCurve(fit1), conf.int = TRUE)

## optional: get MCMC samples for more estimates of interest
mcmc_list <- read.bayesCox(fit1)
str(mcmc_list)

### Fit dynamic coefficient model with time-varying baseline hazards
fit2 <- bayesCox(myformula, mydata, out = "dynCox1.txt", model = "Dynamic",
                 base.prior = list(type = "Gamma", shape = 0.1, rate = 0.1),
                 coef.prior = list(type = "HAR1", shape = 2, scale = 1),
                 gibbs = list(iter = 100, burn = 20, thin = 1,
                              verbose = TRUE, nReport = 20))
plotCoef(coef(fit2))
plotJumpTrace(jump(fit2))
plotJumpHist(jump(fit2))
plotNu(nu(fit2))
plotSurv(survCurve(fit2), conf.int = TRUE)

## Plot the coefficient estimates from three models together
plotCoef(rbind(coef(fit0), coef(fit1), coef(fit2)))

## optional: get MCMC samples for more estimates of interest
mcmc_list <- read.bayesCox(fit2)
str(mcmc_list)

### Fit dynamic coefficient model with dynamic hazards (in log scales)
fit3 <- bayesCox(myformula, mydata, out = "dynCox2.txt", model = "Dynamic",
                 base.prior = list(type = "Const"),
                 coef.prior = list(type = "HAR1", shape = 2, scale = 1),
                 gibbs = list(iter = 100, burn = 20, thin = 1,
                              verbose = TRUE, nReport = 20),
                 control = list(intercept = TRUE))
plotCoef(coef(fit3))
plotJumpTrace(jump(fit3))
plotJumpHist(jump(fit3))
plotNu(nu(fit3))
plotSurv(survCurve(fit3), conf.int = TRUE)

## Plot the estimated survival function and the difference
plotSurv(survCurve(fit3, newdata = newDat, type = "survival"),
         legendName = "Treatment", conf.int = TRUE)
plotSurv(survDiff(fit3, newdata = newDat, type = "survival"),
         legendName = "Treatment", conf.int = TRUE, smooth = TRUE)

## optional: get MCMC samples for more estimates of interest
mcmc_list <- read.bayesCox(fit3, parts = "all")
str(mcmc_list)

}

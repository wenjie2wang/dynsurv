\dontrun{

library(dynsurv)
set.seed(1216)

## breast cancer data
data(bcos)
mydata <- bcos
myformula <- Surv(left, right, type = "interval2") ~ trt

### Fit time-independent coefficient model
fit0 <- bayesCox(myformula, mydata, model = "TimeIndep",
                 base.prior = list(type = "Gamma", shape = 0.1, rate = 0.1),
                 coef.prior = list(type = "Normal", mean = 0, sd = 1),
                 gibbs = list(iter = 100, burn = 20, thin = 1, verbose = FALSE))

## plot coefficient estimates
plotCoef(coef(fit0, level = 0.9))

## Plot the estimated survival function for given new data
newDat <- data.frame(trt = c("Rad", "RadChem"))
row.names(newDat) <- unique(newDat$trt)
plotSurv(survCurve(fit0, newDat), conf.int = TRUE)

### Fit time-varying coefficient model
fit1 <- bayesCox(myformula, mydata, model = "TimeVary",
                 base.prior = list(type = "Gamma", shape = 0.1, rate = 0.1),
                 coef.prior = list(type = "AR1", sd = 1),
                 gibbs = list(iter = 100, burn = 20, thin = 1,
                              verbose = TRUE, nReport = 20))

plotCoef(coef(fit1))
plotSurv(survCurve(fit1), conf.int = TRUE)

### Fit dynamic coefficient model with time-varying baseline hazards
fit2 <- bayesCox(myformula, mydata, model = "Dynamic",
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

### Fit dynamic coefficient model with dynamic hazards (in log scales)
fit3 <- bayesCox(myformula, mydata, model = "Dynamic",
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

## extract MCMC samples
mcmc_list <- bayesCoxMcmc(fit3, part = "coef")
posterior_coef <- mcmc_list$coef
## posterior probabilities of hazard ratio of RadChem (vs. Rad)
## greater than 1 at time 10
posterior_coef[covariate == "trtRadChem" & time == 10, mean(exp(coef) > 1)]

}

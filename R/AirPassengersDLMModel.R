# File:   
#
# 
# 
# 
# Inspired on: https://lbelzile.github.io/timeseRies/state-space-models-and-the-kalman-filter.html
# INSTALL AND LOAD PACKAGES ################################

# Installs pacman ("package manager") if needed
if (!require("pacman")) install.packages("pacman")

# Use pacman to load add-on packages as desired
# Packages I load every time; uses "pacman"
pacman::p_load(pacman,tidyverse,lubridate,gridExtra,dlm) 



# Current worklog ##########################################











library(datasets)
library(dlm)
data("AirPassengers")
plot(log(AirPassengers))
y <- as.numeric(log(window(AirPassengers,end=c(1958))))
time <- as.vector(time(window(AirPassengers,end=c(1958))))
plot(time,y,type="l")
# Input values came from a MLE and I skip this with estimated values already
dlmModel <- dlmModPoly(order=2,dV=0.0001257008,dW=c(6.291096e-04,8.283651e-10))+
  dlmModSeas(frequency=12,dV=0.0001257008,dW=c(2.20787e-05, rep(0, 10)))
# Variance matrices
V(dlmModel)#; #W(dlmModel), but only diag(W(dlmModel))[1:3] non-null
# Kalman filter
filtered <- dlmFilter(y,dlmModel)
forecasted <- dlmForecast(filtered,nAhead=36)

# list of variances of future observations 
# can be used to create confidence intervals
# data.frame(f=forecasted$f) %>% ggplot()+geom_density(aes(f))
standard_deviations <- sqrt(unlist(forecasted$Q))
# qnorm(0.95) calculates the 95th percentile of the standard normal distribution
# Confidence intervals
plot(log(window(AirPassengers,start=1956)),ylim=c(5,7))
timelo <- seq(tail(time,1)+1/12,by=1/12,length=36)
polygon(x=c(timelo, rev(timelo)),
        y=c(forecasted$f+qnorm(0.95)*standard_deviations,
            rev(forecasted$f-qnorm(0.95)*standard_deviations)),
        lty=3,col=scales::alpha("blue",alpha=0.1))
lines(timelo,forecasted$f,col=3)
legend("bottomright", c("Original series","36-step ahead forecast","Confidence interval"), 
       lty=c(1, 1, 3), lwd=c(2,2,2), col=c(1,3,4),bty = "n")

residuals(filtered)







# Air Passengers dataset ###################################

library(dlm)
data("AirPassengers")
#?"AirPassengers"
#Fit a local level with seasonal dummies
plot(log(AirPassengers), ylab = "log monthly total passengers (in thousands)", 
     main = "Number of international air passengers")
#Keep some observations for forecast validation
y <- as.numeric(log(window(AirPassengers, end = c(1958))))
n <- length(y)
time <- as.vector(time(window(AirPassengers, end = c(1958))))
month <- as.factor(as.integer(time*12) %%12)
#Build the dynamic linear model  - specification with local level + seasonal component (dummies)
#Need to specify initial values for mean m0, otherwise forecasts are diffuse values from zero
#takes a few steps to adjust
build <- function(params) {
  level <- dlmModPoly(order = 2, dV = exp(params[1]), dW = c(exp(params[2:3])))
  seas <- dlmModSeas(frequency = 12, dV = exp(params[1]), dW = c(exp(params[4]), rep(0, 10)))  #stochastic seas.
  mod <- level + seas 
  return(mod)
}
#Initial parameter values - four log-variance, mean and trend, seasonal effect
init <- rep(-5,4)
#Fit the DLM model - numerical optimization
fit <- dlmMLE(y = y, parm = init, build = build, method = "BFGS", control=list(trace = 10, maxit = 1000)) 





#Define model with estimated parameters
mod <- build(fit$par)

#Variance matrices
V(mod)#; #W(mod), but only diag(W(mod))[1:3] non-null







#Filtering
filtered <- dlmFilter(y, mod) 
#plot(time, y, type = "l")
lines(time, c(mod$FF %*% t(filtered$m[-1,])), col = 2, lwd = 2) #filtered states

#One-step ahead forecasts (a linear fn of filtering mean)
forecasted <- dlmForecast(filtered, nAhead = 36)
timelo <- seq(tail(time,1) + 1/12, by = 1/12, length = 36)
lines(timelo,  forecasted$f, col = 4, lwd=2)
legend("bottomright", c("Original series","Filtered states","One-step forecasts"), 
       lty=c(1, 1, 1), lwd=c(1,2,2), col=c(1,2,4),bty = "n")









#90% Confidence intervals
plot(log(window(AirPassengers, start=1956)), ylab = "log monthly total passengers (in thousands)", 
     main = "Number of international air passengers", ylim=c(5,7), type = "b", pch = 20)
polygon(x=c(timelo, rev(timelo)), y = c(forecasted$f + qnorm(0.95)*sqrt(unlist(forecasted$Q)), 
                                        rev(forecasted$f - qnorm(0.95)*sqrt(unlist(forecasted$Q)))), col=scales::alpha("blue", alpha=0.3))
lines(timelo, forecasted$f, col = 4, lty = 2)
legend("topleft", c("Data","One-step forecasts"), lty=c(1, 2), lwd=2, col=c(1,4), pch=c(20,NA), bty = "n")




# Predictions ahead
pred_air <- predict(model_air, interval = "prediction", level = 0.9, n.ahead = 36)
ts.plot(cbind(log(window(AirPassengers, start = 1955, end = 1961)), pred_air), 
        col = c(1, 2, "grey", "grey"), ylab = "log monthly total passengers (in thousands)", 
        main = "Number of international air passengers", lty = c(1, 1, 2, 2))
legend("topleft", c("Observed", "Predictions", "Prediction intervals"), lty = c(1, 
                                                                                1, 2), lwd = 2, col = c(1, 2, "grey"), bty = "n")




smoothed <- dlmSmooth(filtered)
plot(time, smoothed$s[-1,1], type= "l", ylab = "", ylim=c(4,6.5), lwd=2, main = "Smoothed states")
lines(time, 4.5 + smoothed$s[-1,3], col = 2, lwd = 2, lty = 4)
legend("topleft", c("Local linear","Seasonal + 4.5"), lty = c(1, 4), lwd = 2, col = 1:2, bty = "n")
abline(lm(y ~ time)$coefficients)







par(oldpar)

# Normal quantile-quantile plot with confidence intervals
res <- resid(filtered)$res
n <- length(res)
quant <- qnorm(rank(res)/(n + 1))
# Pointwise limits based on distribution of order statistics
confint_lim <- t(sapply(1:n, function(i) {
  qnorm(qbeta(c(0.025, 0.975), i, n - i + 1))
}))
matplot(sort(quant), confint_lim, type = "l", lty = 2, col = "grey", main = "Normal QQ-plot", 
        xlab = "Theoretical quantiles", ylab = "Empirical quantiles")
# Pointwise limits based on draws from null model
confint_lim2 <- t(apply(apply(matrix(rnorm(1000 * n), nrow = n, ncol = 1000), 
                              2, sort), 1, function(x) {
                                quantile(x, c(0.025, 0.975))
                              }))
matplot(sort(quant), confint_lim, type = "l", lty = 2, col = "grey", main = "Normal QQ-plot", 
        xlab = "Theoretical quantiles", ylab = "Empirical quantiles")
matplot(sort(quant), confint_lim2, type = "l", lty = 3, col = "lightgrey", add = TRUE)
points(quant, res, pch = 20)
abline(c(0, 1))







# Ljung-Box test on residuals
plot(sapply(4:24, function(l) {
  Box.test(x = res, lag = l, type = "L", fitdf = 3)$p.value
}), type = "h", ylim = c(0, 1), ylab = "P-value", main = "Ljung-Box test on residuals", 
bty = "l")
abline(h = 0.05, col = 2)








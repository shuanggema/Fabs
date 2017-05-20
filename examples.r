library(lars)
library(glmnet)
source("standardize.r")
source("loss.r")
source("penalty.r")
source("derivative.r")
source("fabs.r")
source("fabsBridge.r")

## For Windows
source("spr.r")
## For Mac,Linux, compile "spr.c" file and then
## replace 'dyn.load("spr.dll")' to 'dyn.load("spr.so")' in "spr.r" file.

# -------------------------------------------------------------------------------------
# source data from packages "lars"
data(diabetes)
 
unAsIs <- function(X) {
  if("AsIs" %in% class(X)) {
    class(X) <- class(X)[-match("AsIs", class(X))]
  }
  X
}

x = unAsIs(diabetes$x)
y = diabetes$y
  
n = nrow(x)
p = ncol(x)
group = 1:p
  
# initialize parameters.
model      = "lm"        # loss function used. Righte now it could be "lm" or "spr". 
weight     = NULL        # weights for adaptive lasso. Defalt is NULL.
lambda.min = 0           # smallest value for lambda, as a fraction of lambda.max.
status     = rep(1, n)   # censoring indicator. 
sigma      = 1/sqrt(n)   # smoothing parameters in sigmiodal functions.
gamm       = NULL        # tuning in Bridge penalty. "1" implies "Lasso", and "2" implies "Bridge".
back       = TRUE        # whether a backward is taken is Fabs. "FALSE" implies Forward Stagewise. 
stoping    = TRUE        # whether a stoping criteria is used. "TRUE" for Fabs.
eps        = 10^-2       # step size in Fabs.
xi         = 0           # used for experiment. Set it to be 0 or < 10^-6.
iter       = 10^6        # maximal iterations.

   
#############################################################################################################
#############################################################################################################
# fabs for linear model and smoothed partial rank estimation.
fit1 = fabs(x, y, status, sigma, group, "lm",  weight, "gLasso", gamm, back, stoping, eps, xi, iter, lambda.min)
fit2 = fabs(x, y, status, sigma, group, "spr", weight, "gLasso", gamm, back, stoping, eps, xi, iter, lambda.min)    
par(mfrow=c(1,2))	
matplot(colSums(abs(fit1$beta[-1,])), t(fit1$beta[-1,]), lty=1, 
        type="l", col=1:p, lwd=2, xlab = expression(paste("||", symbol(beta), "||",sep="")[1]), 
	    ylab = "Coordinates", main="OLS")
matplot(colSums(abs(fit2$beta)), t(fit2$beta), lty=1, 
        type="l", col=1:p, lwd=2, xlab = expression(paste("||", symbol(beta), "||",sep="")[1]), 
	    ylab = "Coordinates", main="SPR")

#############################################################################################################
#############################################################################################################
# general fabs for Lasso and Ridge regressions.
# compare general fabs with coordinate descent (in glmnet).
eps     = 10^-1
model   ="lm"
# OLS + Lasso 
fit1 = fabsBridge(x, y, status, sigma, group, model, "Bridge", T, 1, TRUE, stoping, eps, xi, iter, lambda.min)
fit2 = glmnet(x, y, alpha = 1, lambda.min.ratio = lambda.min)
# OLS + Ridge
fit3 = fabsBridge(x, y, status, sigma, group, model, "Bridge", F, 2, TRUE, stoping, eps, xi, iter, lambda.min)
fit4 = glmnet(x, y, alpha = 0, lambda.min.ratio = lambda.min)
   
	
par(mfrow=c(1,2))
matplot(colSums(abs(fit1$beta[-1,seq(1, fit1$iter, by=10)])), t(fit1$beta[-1,seq(1, fit1$iter, by=10)]), lty=1, 
        type="l", col=2, lwd=2, xlab = expression(paste("||", symbol(beta), "||",sep="")[1]), 
	    ylab = "Coordinates", main="Lasso")
matplot(colSums(abs(fit2$beta)), t(fit2$beta), lty=3, type="l", col=1, lwd=2,
            pch=19, add = TRUE)
matplot(colSums(abs(fit3$beta[-1,seq(1, fit3$iter, by=100)])), t(fit3$beta[-1,seq(1, fit3$iter, by=100)]), lty=1,
 type="l", col=2, lwd=2, xlab = expression(paste("||", symbol(beta), "||",sep="")[1]), ylab = " ", main="Ridge")
matplot(colSums(abs(fit4$beta)),      t(fit4$beta),      lty=3, type="l", col=1, lwd=2,
            pch=19, add = TRUE)

    


    
 
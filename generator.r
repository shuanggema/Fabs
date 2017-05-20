

# 2015/5/16     in NUFE
require(mvtnorm)
generator = function(n,  p,  d, rho, error, tran, censor.rate=NULL)
{
    sigma = outer(1:p, 1:p, FUN = function(x, y) rho^(abs(x - y)))
    
    x = rmvnorm(n, mean = rep(0,p), sigma = sigma)
    
    #nz = c(runif(floor((d-1)/2), 0.35, 0.4), runif(d-1-floor((d-1)/2), -0.4, -0.35))
    nz = c(rep(1, floor((d-1)/2)), rep(-1, d-1-floor((d-1)/2)))
    b = c(1, nz, rep(0, p-d))
    
    # error          
    if(error == "norm") e = rnorm(n, 0, 1)
    if(error =="contaminate") e = c(rnorm(0.7*n), rcauchy(0.3*n))
    if(substr(error, 1, 1) == "t") e = rt(n, as.numeric(substr(error,2,2)))
    #  the natural log of a Weibull random time is an extreme value random observation.   
    if(error == "ev") e = log(rweibull(n, shape = 1, scale = 1)) - digamma(1) 
    
    # transformation
    if(tran == "x^3")  g = function(x) sign(x) * (abs(x)^(3))
    if(tran == "lm") g = function(x) x
    if(tran == "log") g = function(x) exp(x)
    
    y0 = g(c(x[, 1:d] %*% b[1:d]) + e)
    
    # censor.mu = 0.1, censor rate = 0.16
    #           = 0.2, censor rate = 0.25  
    if   (is.null(censor.rate) != TRUE) {
          cens = quantile(y0, censor.rate)
          y = pmax(y0, cens)
          status = 1 * (y0>=cens)
    }  else {
          y = y0
          status = rep(1, n)
    }
  
    # In the transformation model, y may be infinity.
    y[which(y=="NaN"|y=="Inf")] =  max(y[which(y!="NaN"&y!="Inf")]) 
    list(x=x, y=y, b=b, status=status)                        
}


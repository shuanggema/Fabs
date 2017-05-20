# this is a interface function for compute the  parital rank correlation and
# the smoothed (using sigmoid function ) partial rank correlation.

# sort the data fisrtly, then we can iterate n^2/2 instead of n^2 times
# to calculate the loss functions.
# For Mac/Linux, compile "spr.c" with 
###        system("R CMD SHLIB spr.c") 


dyn.load("spr.dll")

pr <- function(y, x, status, Beta0)
{
    n = length(y)
    p = ncol(x)
    fit <- .C("pr",
              Y = as.double(y),
              X = as.double(c(t(x))),
              Beta = as.double(Beta0),
              correlation = as.double(0),
              Delta = as.integer(status),
              param = as.integer(c(n,p)) )
    list(cor=fit$correlation)
}

hpr <- function(y, x, status, Beta0)
{
    n = length(y)
    p = ncol(x)
    fit <- .C("hpr",
              Y = as.double(y),
              X = as.double(c(t(x))),
              Beta = as.double(Beta0),
              correlation = as.double(0),
              Delta = as.integer(status),
              param = as.integer(c(n,p)) )
    list(cor=fit$correlation)
}


spr <- function(y, x, status, Beta0, sigma)
{
    n = length(y)
    p = ncol(x)
    fit <- .C("spr",
              Y = as.double(y),
              X = as.double(c(t(x))),
              Beta = as.double(Beta0),
              correlation = as.double(0),
              Delta = as.integer(status),
              Sigma = as.double(sigma),
              param = as.integer(c(n,p)) )
    list(cor=fit$correlation)
}

dspr <- function(y, x, status, Beta0, sigma)
{
    n = length(y)          
    p = ncol(x)
    fit <- .C("dspr",
              Y = as.double(y),
              X = as.double(c(t(x))),
              Beta = as.double(Beta0),
              derivative = as.double(rep(0, p)),
              Delta = as.integer(status),
              Sigma = as.double(sigma),
              param = as.integer(c(n,p)) )
    list(d = fit$derivative)
}


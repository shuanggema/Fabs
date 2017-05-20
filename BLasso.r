# algorithm structure
source("difference.r")

BLasso = function(x, y, status=NULL, sigma=NULL, group, model="lm", w, back=TRUE, stoping=TRUE, b0=NULL, eps = 10^-3, xi = 10^-6, iter=10^4, lambda.min={if (nrow(x) >
ncol(x)) 1e-4 else .05})
{
    # standardize x, y
    std = standardize(x)
    XX = std[[1]]
    center = std[[2]]
    scale = std[[3]]
    yy = y  - mean(y) 

    # step1: initialization
    n = nrow(XX)
    p = ncol(XX)
    backward = forward = 0
    
    b = matrix(0, ncol=iter, nrow=p)
    b0 = b[,1]
    lambda = numeric(iter)
    dif = difference(XX, yy, status, model, sigma, b0, eps, 1:p, w)
    d = dif$d
    k = dif$k
    b[,1] = dif$b.new
    #print(dif)
    
    loss0 = Loss(x=XX, y=yy, status=status, b=b0, sigma=sigma, model=model)
    loss = Loss(x=XX, y=yy, status=status,  b=b[,1], sigma=sigma, model=model)
    lambda[1] =  (loss0-loss)/eps
    set = k

    #sigma = quantile(abs(transfer(yy,XX)$X %*% b0)/5, probs=0.01)
    # step2: backward and forward
    for (i in 1:(iter-1))
    {
        # backward
        der1 = difference(XX, yy, status, model, sigma, b[,i], eps, set, w)
        #print(der1)
        d = der1$d[,2] * (sign(b[set,i])==1) +  der1$d[,1] * (sign(b[set,i]) == -1)
        k = set[which.min(d)]
        # force a backward
        b.temp = b[, i]; b.temp[k] = b[k,i] - eps*sign(b[k,i])/w[k]
        loss0 = Loss(x=XX, y=yy, status=status, b=b[,i], sigma=sigma, model=model)
        loss = Loss(x=XX, y=yy, status=status,  b=b.temp, sigma=sigma, model=model)
        
        descent =  loss - loss0 + lambda[i]*(abs(b.temp[k]) - abs(b[k,i]))*w[k]
        if (descent < -xi & back==TRUE) {
            b[,i+1] = b.temp
            lambda[i+1] = lambda[i]
            backward = backward + 1
            if(b[k]==0) set = setdiff(set, k)
        }  else  {
           # forward
           der2 = difference(XX, yy, status, model, sigma, b[,i], eps, 1:p, w) 
           k = der2$k
           b[, i+1] = der2$b.new       ###

           loss0 = Loss(x=XX, y=yy, status=status, b=b[,i], sigma=sigma, model=model)
           loss = Loss(x=XX, y=yy, status=status, b=b[,i+1], sigma=sigma, model=model)
           
           temp = (loss0 - loss - xi)/eps
           
           lambda[i+1] = min(lambda[i], temp)
           #lambda[i+1] = temp
           set = union(set, k)
           #print(set)
           if(b[k, i+1]==0) set = setdiff(set, k)
           forward = forward + 1
        }
        if(stoping==TRUE & lambda[i+1] < lambda[1]*lambda.min) break
        
        if(i==(iter-1)) warning("Algorithm failed to converge.")
        
    }
    print(i)
    
    ## Unstandardize coefficients for the linear model.
    if (model=="lm") 
    {   
      bet = unstandardize(b[, 1:i], center, scale, mean(y)) 
    }  else {
      bet = b[, 1:i]
    }
    
    lambda = lambda[1:i]


    val = list(beta=as.matrix(bet),
               iter=i,
               backward=backward,
               forward=forward,
               lambda=lambda)
    val
}

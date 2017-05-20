# algorithm structure
# 10/6/2015 in NUFE.
# Xingjie Shi
fabs = function(x, y, status=NULL, sigma=NULL, group, model="lm", weight=NULL, type = "gLasso", gamm = NULL,
                back=TRUE, stoping=TRUE, eps = 10^-3, xi = 10^-6,
                iter=10^4, lambda.min={if (nrow(x) > ncol(x)) 1e-4 else .05})
{

    ## Reorder groups, if necessary
    gf    = as.factor(group)
    g     =  as.numeric(gf)
    G     = max(g)
    g.ord = order(g)
    g     = g[g.ord]
    x     = x[,g.ord]
    K     = as.numeric(table(g))
    K1    = cumsum(K)
    K0    = c(1, K1[-G]+1)

    ## set up y, x, w
    std    = standardize(x)
    XX     = std[[1]]
    center = std[[2]]
    scale  = std[[3]]
    yy     = y  - mean(y)
    
    if(is.null(weight)) w = sqrt(K) else w = weight
    
    lambda = direction = numeric(iter)
    
    # step1: initialization (forward step).
    # --------------------------------------------------------------------------
    n = nrow(XX)
    p = ncol(XX)
    b0 = rep(0, p)
    b  = matrix(0, ncol=iter, nrow=p)
    
    der = derivative(XX, yy, status, model, sigma, type, gamm, K0, K1, b0, w, 1:G)
    d   = der$d
    k   = der$k.f             
    b[K0[k]:K1[k], 1] = - eps*dpen(d[K0[k]:K1[k]], type, gamm)/w[k]
    
    loss0 = Loss(XX, yy, status, b0,    sigma, model)
    loss  = Loss(XX, yy, status, b[,1], sigma, model)
    
    rho.dif = w[k] * (pen(b[,1], type, gamm) - pen(b0, type, gamm))
    lambda[1]    =  (loss0-loss)/rho.dif
    direction[1] = 1
    set          = k

    # step2: backward and forward
    # --------------------------------------------------------------------------
    for (i in 1:(iter-1))
    {
        # backward direction.
        der1 = derivative(XX, yy, status, model, sigma, type, gamm, K0, K1, b[,i], w, set)
        d = der1$d
        k = set[which.min(-d[set]*sign(b[set,i]))]
        
        # force a backward
        Delta.k = dpen(d[K0[k]:K1[k]], type, gamm)/w[k]
        b.temp = b[, i]
        b.temp[K0[k]:K1[k]] = b[K0[k]:K1[k], i] + eps*Delta.k
        loss0 = Loss(XX, yy, status, b[,i],  sigma, model)
        loss  = Loss(XX, yy, status, b.temp, sigma, model)
        
        rho.dif = pen(b.temp[K0[k]:K1[k]], type, gamm) - 
                     pen(b[K0[k]:K1[k], i], type, gamm)
        descent =  loss - loss0 + lambda[i] * w[k] * rho.dif
        
        if (round(descent, 6) < -xi & back==TRUE) {
            b[,i+1] = b.temp
            lambda[i+1] = lambda[i]
            direction[i+1] = -1
            if(norm(b[K0[k]:K1[k], i+1], "2") < 1e-6) set = setdiff(set, k)
        }  else  {
           # forward
           der2 = derivative(XX, yy, status, model, sigma, type, gamm, K0, K1, b[,i], w, 1:G)
           d = der2$d
           k = der2$k.f
          
           Delta.k = -dpen(d[K0[k]:K1[k]], type, gamm)/w[k]

           b[, i+1] = b[, i]
           b[K0[k]:K1[k], i+1] = b[K0[k]:K1[k], i] + eps*Delta.k

           loss0 = Loss(XX, yy, status, b[,i],   sigma, model)
           loss  = Loss(XX, yy, status, b[,i+1], sigma, model)
           
           
           rho.dif = w[k]*(pen(b[K0[k]:K1[k], i+1], type, gamm) - 
                         pen(b[K0[k]:K1[k], i],   type, gamm))
           
           temp  = (loss0 - loss - xi)/abs(rho.dif)
           lambda[i+1] = min(lambda[i], temp)

           set = union(set, k)
           if(norm(b[K0[k]:K1[k], i+1], "2") < 1e-6) set = setdiff(set, k)
           direction[i+1] = 1
        }
        if (stoping==TRUE & lambda[i+1] <= lambda[1]*lambda.min) break
        if (i==(iter-1)) 
            warning("Solution path unfinished, more iterations is needed.")
    }
    print(i)
    
    # prepare output.
    # --------------------------------------------------------------------------
    ## Unstandardize coefficients for the linear model.
    if (model=="lm") 
    {   
        bet = unstandardize(b[, 1:i], center, scale, mean(y)) 
    }  else if (model=="spr"){
         bet = b[, 1:i]/scale
    }
    lambda    = lambda[1:i]
    direction = direction[1:i]

    val = list(beta      = bet,
               iter      = i,
               direction = direction,
               group     = group,
               lambda    = lambda)
    val
}

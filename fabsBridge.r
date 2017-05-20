# algorithm structure
# 10/6/2015 in NUFE.
# Xingjie Shi
fabsBridge = function(x, y, status=NULL, sigma=NULL, group, model="lm", type = "gLasso", sparsity=TRUE,
                gamm = NULL, back=TRUE, stoping=TRUE, eps = 10^-3, xi = 10^-6,
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
  
  w = sqrt(K)
  
  lambda = descent = direction = numeric(iter)
  
  # step1: initialization (forward step).
  # --------------------------------------------------------------------------
  n = nrow(XX)
  p = ncol(XX)
  b0 = rep(0, p)
  b  = matrix(0, ncol=iter, nrow=p)
  
  der = derivative(XX, yy, status, model, sigma, type, gamm, K0, K1, b0, w, 1:G)
  d   = der$d
  k   = der$k.f             
  b[K0[k]:K1[k], 1] =  -eps*sign(d[k])/w[k]
  
  loss0 = Loss(XX, yy, status, b0,    sigma, model)
  loss  = Loss(XX, yy, status, b[,1], sigma, model)
  
  rho.dif = pen(b[,1], type, gamm) - pen(b0, type, gamm)
  lambda[1]    =  (loss0-loss)/rho.dif/w[k]
  direction[1] = 1
  descent[1]   = 0
  # differentiable set (also is the nonzero set), used for sparsity penalty.
  if(sparsity == FALSE) set = 1:G else set = k

  # step2: backward and forward
  # --------------------------------------------------------------------------
  for (i in 1:(iter-1))
  {
    # backward direction.
    der1 = derivative(XX, yy, status, model, sigma, type, gamm, K0, K1, b[,i], w, 1:G)
    d = der1$d
    #
    v = d+lambda[i]* dpen(b[,i], type, gamm)
    k = set[which.max(abs(v[set]))]
	Delta.k = -sign(v[k])/w[k]
	
	if ((length(set) != G) & (sparsity == TRUE)) {
	    h = which.max(abs(d[-set]))
		if (-abs(d[h])*eps + lambda[i]*pen(eps, type, gamm) < -abs(v[k])*eps) {
		    k = h
			Delta.k = -sign(d[h])/w[k]
		}
	}
    # force a backward
    b.temp = b[, i]
    b.temp[K0[k]:K1[k]] = b[K0[k]:K1[k], i] + eps*Delta.k
	
    loss0 = Loss(XX, yy, status, b[,i],  sigma, model)
    loss  = Loss(XX, yy, status, b.temp, sigma, model)
    
    rho.dif = pen(b.temp[K0[k]:K1[k]], type, gamm) - 
      pen(b[K0[k]:K1[k], i], type, gamm)
    descent[i] =  loss - loss0 + lambda[i] * w[k] * rho.dif
    if (descent[i] < -xi & back==TRUE) {
      b[,i+1] = b.temp
      lambda[i+1] = lambda[i]
      direction[i+1] = -1

    }  else  {
      # forward
      der2 = derivative(XX, yy, status, model, sigma, type, gamm, K0, K1, b[,i], w, 1:G)
      d = der2$d
      k = der2$k.f

      Delta.k = -sign(d[k])*eps
      b[, i+1] = b[, i]
      b[K0[k]:K1[k], i+1] = b[K0[k]:K1[k], i] + eps*Delta.k
      
      loss0 = Loss(XX, yy, status, b[,i],   sigma, model)
      loss  = Loss(XX, yy, status, b[,i+1], sigma, model)
      
      
      rho.dif = abs(pen(b[K0[k]:K1[k], i+1], type, gamm) - 
                      pen(b[K0[k]:K1[k], i],   type, gamm))
      
      temp  = (loss0 - loss)/rho.dif/w[k]
      lambda[i+1] = min(lambda[i], temp)
      #lambda[i+1] = ifelse(lambda[i]>temp, temp, (lambda[i]+temp)/2)
      if(sparsity == TRUE) set = union(set, k)
      #if(norm(b[K0[k]:K1[k], i+1], "2") == 0) set = setdiff(set, k)
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
    bet = b[, 1:i]
  }
  lambda    = lambda[1:i]
  direction = direction[1:i]
  descent   = descent[1:i]
  val = list(beta      = bet,
             iter      = i,
             direction = direction,
             descent   = descent,
             group     = group,
             lambda    = lambda)
  val
}

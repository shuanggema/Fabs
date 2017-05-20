# tuning selection using BIC 
# 3/30/2016 @ nufe

# imput 
#      y
#      x
#      status
#      b --- solution path 
#      sigma
# output
#      BI  --- bayesian information
#      opt --- optimal tuning

bic = function(x, y, status, b, sigma) 
{
    grid.n = ncol(b)
    n = length(y)
    BI = numeric(length(grid.n))
    for (i in 1:grid.n)
    {
         score = spr(y, x, status, b[, i], sigma)$cor/n/(n-1)
         BI[i] = log(score) - sum(b[,i]!=0)*log(n)/(2*n)
    }
    list(BI  = BI, 
         opt = which.max(BI))

}


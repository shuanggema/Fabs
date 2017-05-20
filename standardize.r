# standardize design matrix.
# from packages {ncvreg}.
standardize = function(x)
{
   p = ncol(x)
   n = nrow(x)
   center = colMeans(x)
   x.mean = x - matrix(rep(center, n), n, p, byrow=T)
   scale = sqrt(colSums(x.mean^2)/n)
   xx = t(t(x.mean)/scale)
   list(xx, center, scale)
}


# unstandardize 
unstandardize = function(b, center, scale, ybar)
{
  beta = matrix(0, nrow=nrow(b)+1, ncol=ncol(b))
  beta[-1,] = b / scale
  beta[1,] = ybar  - c(t(center) %*% beta[-1,])
  beta
}

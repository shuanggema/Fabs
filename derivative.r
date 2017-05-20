# a function to calculate derivatives and find optimal index k in the active set.
# for spr estimator, no anchor is selected.
derivative = function(X, y, status=NULL, model="lm", sigma=NULL, type, gamm, K0, K1, b, w, set)
{
   if(model=="lm") d = - t(X) %*% (y- c(X%*%b))/length(y)
   if(model=="spr") d = -dspr(y, X, status, b, sigma)$d/length(y)/(length(y)-1)

   d.gr = numeric(length(set))
   i = 1
   for(j in set)
   {             
       d.gr[i] = pen(d[K0[j]:K1[j]], type, gamm)/w[j]
       i = i + 1
   }
   ord = order(d.gr, decreasing=TRUE)
   k.f = set[ord[1]]
   k.b = set[ord[length(set)]]
   list(d = d, k.f = k.f, k.b = k.b)
}

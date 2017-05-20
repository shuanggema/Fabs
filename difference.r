# a function to calculate derivatives and find optimal index k in the active set.
difference = function(x, y, status=NULL, model, sigma, b, eps, set, w)
{
   # difference: each coordinate in set change epsilon, how does the loss change.
   d.gr = matrix(0, ncol=2, nrow=length(set))
   i = 1      
   for(j in set)                                    
   {
       b.up = b.down = b
       b.up[j] = b[j] + eps/w[j]
       d.gr[i,1] = Loss(x, y, status, b.up, sigma, model)
       
       b.down[j] = b[j] - eps/w[j]
       d.gr[i,2] = Loss(x, y, status, b.down, sigma, model)
       i = i + 1
   }
    
   ord.up = order(d.gr[,1], decreasing=FALSE)
   ord.down = order(d.gr[,2], decreasing=FALSE)
   
   k.up = set[ord.up[1]]
   k.down = set[ord.down[1]]
   
   k = ifelse(d.gr[which(set==k.up),1] < d.gr[which(set==k.down),2], k.up, k.down)
   
   b.new = b; 
   b.new[k] = ifelse(d.gr[which(set==k.up),1] < d.gr[which(set==k.down),2], b[k] + eps/w[k], b[k] - eps/w[k])
   list(d=d.gr, k=k, b.new=b.new)
}

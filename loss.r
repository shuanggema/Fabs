# loss
Loss = function(x, y, status=NULL, b, sigma=NULL, model="lm")
{
    n = nrow(x) 
    if(model=="lm")  val = sum((y -c(x %*% b))^2)/(2*length(y))
    if(model=="spr") val = -spr(y, x, status, b, sigma)$cor/(n*(n-1)/2)
    if(model=="pr")  val = -pr(y, x, status, b)$cor/(n*(n-1)/2)
    val
}
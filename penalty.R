pen = function(b, type = c("gLasso", "Bridge")[1], gamm = NULL){
  ifelse (type == "gLasso", norm(b, "2"), sum(abs(b)^gamm))
}
  

dpen = function(b, type = c("gLasso", "Bridge")[1], gamm = NULL){
  if(type == "gLasso") b/norm(b, "2") else gamm*abs(b)^(gamm-1)*sign(b)
}

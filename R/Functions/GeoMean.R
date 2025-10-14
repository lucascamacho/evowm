# geometric mean function
geomean = function(vector){
  g = exp(mean(log(vector)))
  return(g)
}
Cont2Gaus <- function(iData,total_iteration=5000,stepsize=0.05)
{
  # Data Continuzed Transformation 
  continuz <- ContTran(iData,total_iteration=total_iteration,stepsize=stepsize)
  # Data Gaussianized transformation
  D <- log(continuz+1)
  Gaus = huge.npn(D)
  return(Gaus)
  
}
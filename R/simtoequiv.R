##### transform simulated graph to perfect casual graph

simtoequiv <- function(edgematrix) {
  DataP = dim(edgematrix)[1]
  moralM = edgematrix + t(edgematrix)  
  for(i in 1:DataP) {
    for(j in 1:DataP) {
      if (i != j & moralM[i, j] == 0 ){
        kindx <- intersect(which((edgematrix)[i, ] == 1), which((edgematrix)[j, ] == 1))
        
        if(length(kindx) > 0) {
          for(kk in kindx) {
            #cat("i, j, k", i, j, kk, "\n")
            moralM[kk, i] = 0
            moralM[kk, j] = 0
          }
        }
      }
    }
  }
  return(list(PDAG = moralM))
  
}
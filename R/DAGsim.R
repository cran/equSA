
  
DAGsim <- function(n, p, sparsity = 0.02,  p.binary)  {
  edgematrix <- matrix(0, p, p)
  coef = 0.5
  #### adjacency matrix
  for(i in 1:p){  
    for(j in 1:p){
      if(i  > j) {
        edgematrix[i, j] <- rbinom(1, 1, sparsity)
      }
    }
  }
  
  
  X  = matrix(NA, n, p)
  
  
  #### half binary and half Gaussian
  binary_indx <- sample(c(1:p), p.binary)
  cont_indx <- setdiff(c(1:p), binary_indx)
  
  X <- matrix(0, n, p)
  if(1 %in% cont_indx) {
    X[, 1] = rnorm(n, 0, 1)
  } else {
    temp = rnorm(n, 0, 1)
    X[,1] = rbinom (n, 1, prob =  exp(temp)/ (1 + exp(temp)))
  }
  
  for(i in 2:p) {

      if(i  %in% cont_indx) {
        temp = 0
        for(k in 1:(i-1)) {
          temp =  temp + coef*edgematrix[i, k]*X[, k]
        }
        
        X[, i] = temp + rnorm(n, 0, 1)
        
      } else {
        temp = 0
        for(k in 1:(i-1)) {
          temp =  temp + coef*edgematrix[i, k]*X[, k]
        }
        
        X[, i] = rbinom (n, 1, prob =  exp(temp)/ (1 + exp(temp)))
        
      }
  }
  return(list(Adjacency.matrix = t(edgematrix), data = X, gaussian.index = cont_indx,
              binary.index = setdiff(c(1:p), cont_indx)))
}

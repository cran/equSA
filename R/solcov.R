
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'


### calculate the covariance matrix with known structure
### return an array of 2 matrices
### result[1,,] is cov and resutl[2,,] is the precision matrix
### data is cov(data), which is a p*p matrix

### calculate the distance of matrix in L1
disMat = function(A, B){
  temp = sum(abs(A-B))
  return(temp)
}

### move the last col and row of matrix to up left
updateMat = function(A, dim){
  temp = A[c(dim, 1:dim-1), c(dim, 1:dim-1)]
  return(temp)
}

### calculate the covariance matrix with known structure
### return an array of 2 matrices
### result[1,,] is cov and resutl[2,,] is the precision matrix
### data is cov(data), which is a p*p matrix

solcov <-  function(data, struct, tol=10^-5){
  #S = var(data)
  S = cov(data)
  W = S
  p = nrow(S)
  precision = matrix(0,p,p)
  
  while(TRUE){
    temp = W
    
    for(j in 1:p){
      W11 = W[-p,-p]
      structure12 = struct[-p,p]
      S12 = S[-p,p]
      
      ## temp2 record the number of row with connectivity
      if(any(c(1:(p-1))[(structure12==1)])){
        temp2 = c(1:(p-1))[(structure12==1)]
        
        beta1 = solve(W11[temp2,temp2], S12[temp2])
        beta2 = rep(0,p-1)
        beta2[temp2] = beta1
        
        ## update W12
        W[-p,p] = W11 %*% beta2
        W[p,-p] = W[-p,p]
      }
      
      ## update W, structure and S matrix
      W = updateMat(W, p)
      S = updateMat(S, p)
      struct = updateMat(struct, p)
    }
    
    if(disMat(W,temp)<=tol) {
      ## calculate the precision matrix in the final cycle
      for(j in 1:p){
        W11 = W[-p,-p]
        structure12 = struct[-p,p]
        S12 = S[-p,p]
        
        ## temp2 record the number of row with connectivity
        if(any(c(1:(p-1))[(structure12==1)])){
          temp2 = c(1:(p-1))[(structure12==1)]
          
          beta1 = solve(W11[temp2,temp2], S12[temp2])
          beta2 = rep(0,p-1)
          beta2[temp2] = beta1
          
          ## update W12
          W[-p,p] = W11 %*% beta2
          W[p,-p] = W[-p,p]
          
          ## get the last col of precision
          precision[p,p] = 1/(S[p,p] - W[-p,p] %*% beta2)
          precision[-p,p] = -beta2*precision[p,p]
          precision[p,-p] = precision[-p,p]
        }
        
        ## update precision, W, structure and S matrix
        W = updateMat(W, p)
        S = updateMat(S, p)
        struct = updateMat(struct, p)
        precision = updateMat(precision, p)
      }
      
      break
    }
  }
  
  result <- list()
  
  
   result$COV = W
   result$PRE = precision

  
  ## return an array of 2 matrices
  return(result)
}  # end solcov


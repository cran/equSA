

DAGsim <- function(n, p, sparsity = 0.02,  p.binary=round(p/2),type="AR(2)", verbose = TRUE)  {
  if(verbose) cat("Generating data and moral adjacency matrix from", type,"DAG structure....\n")
  
  if(type=="AR(2)"){
    coef = 0.5
    edgematrix =  A_moral = matrix(0, p, p)
    
    #### adjacency matrix
    for(i in 1:p)
    {
      for(j in 1:p)
      {
        if(i<j){
          if(abs(j-i)==1) {edgematrix[i,j]=1;}
          else if(abs(j-i)==2) {edgematrix[i,j]=1;}
        }
      }
    }
    
    
    #### half binary and half Gaussian
    binary_indx <- sample(c(1:p), p.binary)
    cont_indx <- setdiff(c(1:p), binary_indx)
    
    
    ## construct moral graph ##
    e = empty.graph(as.character(1:p))
    amat(e) <- edgematrix
    A_moral <- amat(moral(e))
    
    X <- matrix(0, n, p)
    if(1 %in% cont_indx){
      X[, 1] = rnorm(n, 0, 1)
    } else {
      temp = rnorm(n, 0, 1)
      X[,1] = rbinom (n, 1, prob =  exp(temp)/ (1 + exp(temp)))
    }
    
    for(i in 2:p) {
      
      if(i  %in% cont_indx) {
        temp = 0
        for(k in 1:(i-1)) {
          temp =  temp + coef*edgematrix[k, i]*X[, k]
        }
        
        X[, i] = temp+rnorm(n,0,1)
        
      } else {
        temp = 0
        for(k in 1:(i-1)) {
          temp =  temp + coef*edgematrix[k, i]*X[, k]
        }
        
        X[, i] = rbinom (n, 1, prob =  exp(temp)/ (1 + exp(temp)))
        
      }
    }
  }else if(type=="random"){
    
    coef = 0.5
    edgematrix =  A_moral = matrix(0, p, p)
    
    #### adjacency matrix
    for(i in 1:p){  
      for(j in 1:p){
        if(i  < j) {
          edgematrix[i, j] <- rbinom(1, 1, sparsity)
        }
      }
    }
    
    
    
    
    
    binary_indx <- sample(c(1:p), p.binary)
    cont_indx <- setdiff(c(1:p), binary_indx)
    
    ## construct moral graph ##
    e = empty.graph(as.character(1:p))
    amat(e) <- edgematrix
    A_moral <- amat(moral(e))
    
    X <- matrix(0, n, p)
    if(1 %in% cont_indx){
      X[, 1] = rnorm(n, 0, 1)
    } else {
      temp = rnorm(n, 0, 1)
      X[,1] = rbinom (n, 1, prob =  exp(temp)/ (1 + exp(temp)))
    }
    
    for(i in 2:p) {
      
      if(i  %in% cont_indx) {
        temp = 0
        for(k in 1:(i-1)) {
          temp =  temp + coef*edgematrix[k, i]*X[, k]
        }
        
        X[, i] = temp+rnorm(n,0,1)
        
      } else {
        temp = 0
        for(k in 1:(i-1)) {
          temp =  temp + coef*edgematrix[k, i]*X[, k]
        }
        
        X[, i] = rbinom (n, 1, prob =  exp(temp)/ (1 + exp(temp)))
        
      }
    }
  }else if(type=="alarm"){
    edgematrix <- bnlearn::amat(alarm)
    coef <- 0.5
    
    p <- dim(edgematrix)[1]
    p.binary=round(p/2)
    binary_indx <- sample(c(1:p), p.binary)
    cont_indx <- setdiff(c(1:p), binary_indx)
    
    ## construct moral graph ##
    A_moral <- bnlearn::amat(bnlearn::moral(alarm))
    init_node <- as.numeric(which(apply(edgematrix,2,sum)==0))
    N_node <- length(init_node)
    X <- matrix(0, n, p)
    X[,init_node] <- rmvnorm(n, rep(0,N_node), diag(1,N_node))
    for(i in 1:N_node){
      if(!(init_node[i] %in% cont_indx)) {
        X[, init_node[i]] = rbinom (n, 1, prob =  exp(X[, init_node[i]])/ (1 + exp(X[, init_node[i]])))
      }
    }
    iter_gen <- 1
    while(N_node < p){
      next_node <- setdiff(unique(as.numeric(which(edgematrix[init_node,]==1,arr.ind = TRUE)[,2])),init_node)
      cat("iter=",iter_gen,"\n")
      #cat(next_node,"\n")
      for(i in next_node){
        parent <- as.numeric(which(edgematrix[,i]==1))
        if(all(parent %in% init_node)){
          if(i  %in% cont_indx) {
            temp = 0
            for(k in parent) {
              temp =  temp + coef*edgematrix[k, i]*X[, k]
            }
            
            X[, i] = temp + rnorm(n, 0, 1)
            
          } else {
            temp = 0
            for(k in parent) {
              temp =  temp + coef*edgematrix[k, i]*X[, k]
            }
            
            X[, i] = rbinom (n, 1, prob =  exp(temp)/ (1 + exp(temp)))
            
          }
          init_node <- unique(c(init_node, i))
          N_node <- length(init_node)
        }
      }
      iter_gen = iter_gen +1
    }
    
  }else{
    cat("only 'AR(2)', 'random' and 'alarm' types are provided. \n")
  }
  
  return(list(edgematrix = edgematrix, data = X, moral.matrix = A_moral, gaussian.index = cont_indx,
              binary.index = setdiff(c(1:p), cont_indx)))
}

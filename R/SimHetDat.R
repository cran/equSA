
SimHetDat <- function(n = 100, p = 200, M = 3, mu = 0.3, type = "band") {
  if (type != "band") {
    stop("Only 'band' type is provided !")
  }
  C <- diag(1,p)
  for(i in 1:p)
  {
    for(j in 1:p)
    {
      if(abs(j-i)==1) C[i,j]=0.5
      else if(abs(j-i)==2) C[i,j]=0.25
    }
  }
  
if(M==2){
  mu1_true <- rep(0,p)
  A <- diag(0,p)
  
  sigma <- solve(C)
  
  sigma_true=sigma;
  for(i in 1:p)
  {
    for(j in 1:p)
    {
      if(sigma[i,j]<=sigma[j,i]) {sigma_true[i,j]=sigma[j,i]}
      else {sigma_true[i,j]=sigma[i,j]};
      if(C[i,j]!=0&&i!=j) A[i,j]=1;
    }
  }
  x1 <- rmvnorm(n,mu1_true,sigma_true)
  #### second cluster #####
  mu2_true <- rep(mu,p)
  x2 <- rmvnorm(n,mu2_true,sigma_true)
  data <- rbind(x1,x2)
  label <- c(rep(1,n),rep(2,n))
}else if(M==3){
  mu1_true <- rep(0,p)
  A <- diag(0,p)
  
  sigma <- solve(C)
  
  sigma_true=sigma;
  for(i in 1:p)
  {
    for(j in 1:p)
    {
      if(sigma[i,j]<=sigma[j,i]) {sigma_true[i,j]=sigma[j,i]}
      else {sigma_true[i,j]=sigma[i,j]};
      if(C[i,j]!=0&&i!=j) A[i,j]=1;
    }
  }
  x1 <- rmvnorm(n,mu1_true,sigma_true)
  #### second cluster #####
  mu2_true <- rep(mu,p)
  x2 <- rmvnorm(n,mu2_true,sigma_true)
  mu3_true <- rep(-mu,p)
  x3 <- rmvnorm(n,mu3_true,sigma_true)
  data <- rbind(x1,x2,x3)
  label <- c(rep(1,n),rep(2,n),rep(3,n))
}else{
  stop("Only M=2 or M=3 are provided !")
}
  result <- list(data = data, A = A, label = label)
  return(result)
}

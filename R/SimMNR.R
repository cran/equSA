SimMNR <- function(n,p,coef,family="gaussian"){
  mu_x_true <- rep(0,p)
  sigma_eps_true <- 1
  beta_true <- coef[-1]
  beta_0 <- rep(coef[1],n)
  #beta_true[1:5] <- c(10,20,-15,-25,50)/5
  C <- diag(1,p)
  for(i in 1:p)
  {
    for(j in 1:p)
    {
      if(abs(j-i)==1) C[i,j]=0.5
      else if(abs(j-i)==2) C[i,j]=0.25
    }
  }
  
  
  A <- diag(0,p)
  
  sigma <- solve(C)
  mu <- rep(0,p)
  
  sigma2=sigma;
  for(i in 1:p)
  {
    for(j in 1:p)
    {
      if(sigma[i,j]<=sigma[j,i]) {sigma2[i,j]=sigma[j,i]}
      else {sigma2[i,j]=sigma[i,j]};
      if(C[i,j]!=0&&i!=j) A[i,j]=1;
    }
  }
  
  if(family=="gaussian"){
    x <- rmvnorm(n,mu_x_true,sigma2)
    eps <- rnorm(n,0,sigma_eps_true)
    y <- beta_0 +x%*%beta_true+eps
  }else if(family=="binomial"){
    sum0=sum1=0;
    y <- x <- NULL;
    while(!(sum0==n/2&&sum1==n/2)){
      x_pot <- rmvnorm(1,mu_x_true,sigma2)
      eps <- rnorm(1,0,sigma_eps_true)
      linpred <- beta_0[1] +x_pot%*%beta_true
      prob <- exp(linpred)/(1+exp(linpred))
      runis <- runif(1,0,1)
      y_pot <- ifelse(runis<prob,1,0)
      if(sum0< n/2 && sum1 < n/2){
        if(y_pot==0){
          sum0=sum0+1
        }else{
          sum1=sum1+1
        }
        y <- c(y,y_pot)
        x <- rbind(x,x_pot)
      }else if(sum0==n/2 && sum1 < n/2){
        if(y_pot==1) {y <- c(y,y_pot);x <- rbind(x,x_pot);sum1=sum1+1}
      }else if(sum1==n/2 && sum0 < n/2){
        if(y_pot==0) {y <- c(y,y_pot);x <- rbind(x,x_pot);sum0=sum0+1}
      }
    }
  }else if(family=="cox"){
    x <- rmvnorm(n,mu,sigma2)
    
    lambdaT = 0.1 # baseline hazard
    lambdaC = 1 # hazard of censoring
    
    # true event time
    Time = rweibull(n, shape=1, scale=lambdaT*exp(-x%*%beta_true))
    C = rweibull(n, shape=1, scale=lambdaC)   #censoring time
    y0 = pmin(Time,C)  #observed time is min of censored and true
    delta = as.numeric(y0==Time)   # set to 1 if event is observed
    
    ## MNR method ##
    y <- survival::Surv(y0,delta)
  }else{
    stop("Only 'gaussian', 'binomial' and 'cox' models are provided.\n")
  }
  
  return(list(x=x,y=y,Adj=A))
  

  
}
  
  
  
GGMM <- function(data, A, M, alpha1 = 0.1, alpha2 = 0.05, alpha3 = 0.05, iteration = 30, warm = 20)
{
  leng <- iteration - warm
  n <- nrow(data)
  p <- ncol(data)
  GraRes <- equSAR(data, ALPHA1 = alpha1, ALPHA2 = alpha2)
  Ak <- GraRes$Adj
  ## initial ###
  tau_initial <- NULL;data_all <- NULL;
  ind <- sample(1:M,n,replace=TRUE)
  mu_initial <- matrix(1:(M*p),ncol=M)
  for(k in 1:M)
  {
    tau_initial[k] <- length(which(ind==k))/n
    datax <- data[which(ind==k),]
    mu_initial[,k] <- apply(datax,2,mean)
    par <- dim(datax)
    n_k <- par[1]
    for(j in 1:p)
    {
      datax[,j] <- datax[,j]-rep(mu_initial[,k][j],n_k)
    }
    data_all <- rbind(data_all,datax)
  }
  sol <- solcov(data_all,Ak)
  sigma_initial <- sol$COV
  pre_initial <- sol$PRE
  mu <- mu_initial
  sigma <- sigma_initial
  pre <- pre_initial
  tau <- tau_initial
  
  
  ######### IC #####
  U_last <- NULL;
  y_last <- NULL;tau_hat<- NULL;res_hat <- NULL;
  for(iter in 1:iteration)
  {
    Cfy <- log_Cfy <- matrix(1:(n*M),nrow=n)
    for(i in 1:n)
    {
      dd <- NULL;
      for(k in 1:M)
      {
        dd[k] <- dmvnorm(data[i,],mu[,k],sigma,log=TRUE)
      }
      Msum <- 0;
      for(l in 1:M)
      {
        Msum <- Msum +exp(log(tau[l])+dmvnorm(data[i,],mu[,l],sigma,log=TRUE)-max(dd))
      }
      for(k in 1:M)
      {
        log_Cfy[i,k] <- log(tau[k])+dmvnorm(data[i,],mu[,k],sigma,log=TRUE) - max(dd)- log(Msum)
        Cfy[i,k] <- exp(log_Cfy[i,k])
      }
    }
    y <- NULL;
    for(i in 1:n)
    {
      y[i] <- sample(x=1:M,size=1,replace=TRUE,prob=as.vector(Cfy[i,]))
    }
    U <- NULL; data_all <- NULL;
    #### calculate psi scores ####
    for(k in 1:M)
    {
      tau[k] <- length(which(y==k))/n
      datax <- data[y==k,]
      par <- dim(datax)
      n_k <- par[1]
      mu[,k] <- apply(datax,2,mean) 
      for(j in 1:p)
      {
        datax[,j] <- datax[,j]-rep(mu[,k][j],n_k)
      }
      data_all <- rbind(data_all,datax)
      U1 <- psical(datax, ALPHA1 = alpha1)
      U <- cbind(U, U1[, 3])
    }
    #### combine M cluster ######
    N<-p*(p-1)/2
    ratio<-ceiling(N/100000)
    index <- U1[,1:2]
    U<-cbind(index,U)
    z <- apply(tau*U[,-c(1,2),drop=FALSE],1,sum)/sqrt(sum(tau^2))
    if(iter > warm){
      U_last <- cbind(U_last,z)
      y_last <- rbind(y_last,y)
    }
    q<-pnorm(-abs(z), log.p=TRUE)
    q<-q+log(2.0)
    s<-qnorm(q,log.p=TRUE)
    s<-(-1)*s
    U<-cbind(index,s)
    U[U==-Inf] <- -3+rnorm(1,0,0.1)
    U[U== Inf] <- max(U[,3]) + rnorm(1,2,0.1)
    psiscore <- U
    U<-U[order(U[,3]), 1:3]
    m<-floor(N/ratio)
    m0<-N-m*ratio
    s<-sample.int(ratio,m,replace=TRUE)
    for(i in 1:length(s)) s[i]<-s[i]+(i-1)*ratio
    if(m0>0){
      s0<-sample.int(m0,1)+length(s)*ratio
      s<-c(s,s0)
    }
    Us<-U[s,]
    
    qqqscore <- pcorselR(Us, ALPHA2 = alpha3)
    U <- psiscore
    U <- U[U[, 3] >= qqqscore, ]
    Ak <- matrix(rep(0, p * p), ncol = p)
    for (i in 1:nrow(U))
    {
      k1 <- U[i, 1]
      k2 <- U[i, 2]
      Ak[k1, k2] <- 1
      Ak[k2, k1] <- 1
    }
    
    sol <- solcov(data_all,Ak)
    sigma<- sol$COV
    pre <- sol$PRE
    
  }

  
  ######  use last U_last to get A^f ######
  
  U_last<-cbind(index,U_last)
  z <- apply(as.matrix(U_last[,-c(1,2),drop = FALSE]),1,sum)/leng
  
  q<-pnorm(-abs(z), log.p=TRUE)
  q<-q+log(2.0)
  s<-qnorm(q,log.p=TRUE)
  s<-(-1)*s
  U_last<-cbind(U_last[,1:2],s)
  U_last[U_last==-Inf] <- -3+rnorm(1,0,0.1)
  U_last[U_last== Inf] <- max(U_last[,3]) + rnorm(1,2,0.1)
  
  ### calculate PreRec ###
  fit <- summary(U_last[, 3])
  cut1 <- seq(fit[1], fit[2], length = 10)
  cut2 <- seq(fit[2], fit[3], length = 10)
  cut3 <- seq(fit[3], fit[5], length = 10)
  cut4 <- seq(fit[5], fit[6], length = 30)
  cut <- c(cut1, cut2, cut3, cut4)
  precision <- recall <- NULL
  w3 <- NULL
  for (i in 1:length(cut))
  {
    A2 <- diag(0, p)
    aa <- U_last[cut[i] - U_last[, 3] < 0.001, ,drop = FALSE]
    for (j in 1:nrow(aa))
    {
      A2[aa[j, 1], aa[j, 2]] <- 1
      A2[aa[j, 2], aa[j, 1]] <- 1
    }
    B <- A + A2
    L <- length(which(B == 2))
    N <- length(which(B == 0))
    precision <- 1 * L/sum(A2)
    recall <- 1 * L/sum(A)
    w3 <- rbind(w3, c(recall, precision))
  }
  RecPre <- w3
  colnames(RecPre) <- c("Recall", "Precision")
  
  psiscore <- U_last
  U_last<-U_last[order(U_last[,3]), 1:3]
  N<-p*(p-1)/2
  ratio<-ceiling(N/100000)
  m<-floor(N/ratio)
  m0<-N-m*ratio
  s<-sample.int(ratio,m,replace=TRUE)
  for(i in 1:length(s)) s[i]<-s[i]+(i-1)*ratio
  if(m0>0){
    s0<-sample.int(m0,1)+length(s)*ratio
    s<-c(s,s0)
  }
  Us<-U_last[s,]
  
  qqqscore <- pcorselR(Us, ALPHA2 = alpha3)
  
  U <- psiscore
  U <- U[U[, 3] >= qqqscore, ]
  Adj <- matrix(rep(0, p * p), ncol = p)
  for (i in 1:nrow(U))
  {
    k1 <- U[i, 1]
    k2 <- U[i, 2]
    Adj[k1, k2] <- 1
    Adj[k2, k1] <- 1
  }
  

  
### calculate BIC score ##
  loglf_IC <- function(M,mu,sigma)
  {
    Nsum <- 0;
    for(i in 1:n)
    {
      dd <- NULL;
      for(k in 1:M)
      {
        dd[k] <- dmvnorm(data[i,],mu[,k],sigma,log=TRUE)
      }
      Msum <- 0;
      for(l in 1:M)
      {
        Msum <- Msum +exp(log(tau[l])+dmvnorm(data[i,],mu[,l],sigma,log=TRUE)-max(dd))
      }
      Nsum <- Nsum + max(dd)+ log(Msum)
    }
    return(Nsum);
  }
  
  
  matrix_df <- function(w)
  {
    count=0
    for(i in 1:p){
      for(j in i:p)
      {
        if(w[i,j]!=0) count=count+1
      }
    }
    return(count)
  }
  
  
  BIC_IC <- function(A,M,mu,sigma)
  {
    df <- M*(p+ matrix_df(A))
    
    bic <- -2*loglf_IC(M,mu,sigma)+log(n)*df
    return(bic)
  }
  
  BIC_score <- BIC_IC(Adj,M,mu,sigma) 
  
  
  ### identify labels ###
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }

  label <- as.numeric(apply(y_last,2,getmode))
  
  result <- list()
  result$Adj <- Adj
  result$RecPre <- RecPre
  result$BIC <- BIC_score
  result$label <- label
  
  return(result) 
  
}

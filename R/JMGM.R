
JMGM <- function(data,ALPHA1=0.05,ALPHA2=0.01,restrict=FALSE,parallel=FALSE,nCPUs){
  U2 <- NULL;
  M=length(data);
  if(parallel){
    cores <- nCPUs
    cl<-makeCluster(cores)
    registerDoParallel(cl,cores=cores)
    U2 <- foreach(i=1:M,.combine="cbind")  %dopar%
    {
      data_i <- data[[i]]
      score <- plearn.struct(data_i,alpha1=ALPHA1,alpha2=ALPHA2,restrict = restrict, score.only=TRUE)
      if(i==1){
        score
      }else{
        score[,3]
      }
    }
    stopImplicitCluster()
    stopCluster(cl)
  }else{
   for(i in 1:M){
      data_i <- data[[i]]
      score <- plearn.struct(data_i,alpha1=ALPHA1,alpha2=ALPHA2,restrict = restrict, score.only=TRUE)
      if(i==1){
        U2 <- cbind(U2,score)
      }else{
        U2 <- cbind(U2,score[,3])
      }
   }
  }
  score_raw <- U2

### pre_function ####
Ind2 <- function(x){
  e_mod <- ifelse(mean(x)>=0.5,1,0)
  cc <- abs(x-e_mod)
  k1 <- sum(cc)
  k2 <- M-k1
  return(c(k1,k2))
}

#### setting ####
N <- dim(score_raw)[1]
p <- dim(data[[1]])[2]

#### hyperparameter r=2 ###

a1 <- 1
b1 <- 10
a2 <- 1
b2 <- 1

##### S #### 
permu <- function(n){
  S <- NULL;
  for(m in 0:n)
  {
    S <- rbind(S,t(apply(combn(1:n,m=m),2,function(cm) replace(rep(0,n),cm,1))))
  }
  return(S)
}
S <- permu(M)
U <- score_raw
psi_final_hat <- NULL;
num <- dim(S)[1]

for(kk in 1:N)
{
  psi <- as.numeric(U[kk,-c(1,2)])
  wi <- NULL;
  psi_com_hat <- NULL;
  for(ite in 1:num)
  {
    n0 <- length(which(S[ite,]==0))
    n1 <- length(which(S[ite,]==1))
    k1 <- Ind2(S[ite,])[1]  ## change
    k2 <- M-k1  ## not change
    psi0 <- psi[S[ite,]==0]
    psi1 <- psi[S[ite,]==1]
    ## combined psi score ##
    psi_com <- NULL;
    psi0_com <- sum(psi0)/sqrt(length(psi0))
    psi1_com <- sum(psi1)/sqrt(length(psi1))
    psi_com[S[ite,]==0] <- psi0_com
    psi_com[S[ite,]==1] <- psi1_com
    psi_com_hat <- rbind(psi_com_hat,psi_com)
    #### calculate Prob ####
    p0 <- 1/sqrt(n0)*(1/sqrt(2*pi))^n0*gamma(n0/2-1/2+a2)*(1/2*sum(psi0^2)-1/(2*n0)*(sum(psi0))^2+b2)^(-(n0/2-1/2+a2))
    p1 <- 1/sqrt(n1)*(1/sqrt(2*pi))^n1*gamma(n1/2-1/2+a2)*(1/2*sum(psi1^2)-1/(2*n1)*(sum(psi1))^2+b2)^(-(n1/2-1/2+a2))
    if(n0==0)
    {
      wi[ite] <- beta(k1+a1,k2+b1)*p1
    }else if(n1==0){
      wi[ite] <- beta(k1+a1,k2+b1)*p0
    }else{
      wi[ite] <- beta(k1+a1,k2+b1)*p0*p1
    }
  } 
  Prob <- wi/sum(wi) 
  psi_final <- t(t(psi_com_hat)%*%Prob)
  psi_final_hat <- rbind(psi_final_hat,psi_final)
}
score_com <- cbind(score_raw[,c(1,2)],psi_final_hat)


### combine all psi score ###
psi_all <- as.vector(as.matrix(score_com[,-c(1,2)]))
num <- 1:length(psi_all)
U <- cbind(num,num,psi_all)

# subsample
U<-U[order(U[,3]), 1:3]
N<-length(U[,1])
ratio<-ceiling(N/100000)
m<-floor(N/ratio)
m0<-N-m*ratio
s<-sample.int(ratio,m,replace=TRUE)
for(i in 1:length(s)) s[i]<-s[i]+(i-1)*ratio
if(m0>0){
  s0<-sample.int(m0,1)+length(s)*ratio
  s<-c(s,s0)
}
Us<-U[s,]

y <- round(Us,6)
qscore2 <- pcorselR(y,ALPHA2=ALPHA2)


#### cut-off psi score qscore2 ###
A <- array(0,dim = c(M,p,p))
for(k in 1:M)
{
  psik <- score_com[,c(1,2,k+2)]
  psik <-psik[psik[,3]>=qscore2,]
  for(i in 1:nrow(psik)){
    k1<-psik[i,1]
    k2<-psik[i,2]
    A[k,,][k1,k2]<-1
    A[k,,][k2,k1]<-1
  }
}
result <- list()

result$Array = A
result$score.sep = score_raw
result$score.joint = score_com
return(result)

}

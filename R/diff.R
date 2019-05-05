diffR <- function(Data1,Data2,ALPHA1=0.05,ALPHA2=0.05)
  
{

  DataP=length(Data1[1,]);
  
  
  
  ################ Data generation #####################
  p<-DataP;
  N<-p*(p-1)/2;
  ratio<-ceiling(N/100000)
  U1 <- psical(Data1,ALPHA1=ALPHA1);
  U2 <- psical(Data2,ALPHA1=ALPHA1);
  
  U<--(U1[,3]-U2[,3])/sqrt(2)
  U<-cbind(U1[,1:2],U)
  
  z<-U[,3]
  q<-pnorm(-abs(z), log.p=TRUE)
  q<-q+log(2.0)
  s<-qnorm(q,log.p=TRUE)
  s<-(-1)*s
  U<-cbind(U[,1:2],s)
  score <- round(U,6)
  
  ##################### FDR control #################################
  
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
  y <- round(Us,6)
  qqqscore <- pcorselR(y,ALPHA2=ALPHA2);
  
  U<-score
  U<-U[U[,3]>=qqqscore,]
  
  
  A<-matrix(rep(0,p*p), ncol=p)
  for(i in 1:nrow(U)){
    k1<-U[i,1]
    k2<-U[i,2]
    A[k1,k2]<-1
    A[k2,k1]<-1
  }
  return(A)
}

 equSAR <- function(iData,iMaxNei=as.integer(iDataNum/log(iDataNum)),ALPHA1=0.05,ALPHA2=0.05,GRID=2,iteration=100)
{
  iDataNum=length(iData[,1]);
  iDataP=length(iData[1,]);
  Data=NULL;
  for(i in 1:iDataNum)
  {
    for(j in 1:iDataP)
      Data[(i-1)*iDataP+j]=iData[i,j];
  }
  long_mat <- vector("integer",iDataP*iDataP)
  long_score <- vector("double",3*iDataP*(iDataP-1)/2)
   result <- .C("equSA1",
      as.numeric(Data),
      as.integer(iMaxNei),
      as.integer(iDataNum),
      as.integer(iDataP),
      as.numeric(ALPHA1),
      as.numeric(ALPHA2),
      as.integer(GRID),
      as.integer(GRID),
      as.integer(iteration),
      as.integer(long_mat),
      as.numeric(long_score)
      )
   
   A <- matrix(result[[10]],ncol=iDataP,byrow=TRUE)
   score <-matrix(result[[11]], ncol=3, byrow=TRUE)
   result <- list()
   result$Adj = A
   result$score = score
   return(result)
}

psical <- function(iData,iMaxNei=as.integer(iDataNum/log(iDataNum)),ALPHA1=0.05,GRID=2,iteration=100)
{
  iDataNum=length(iData[,1]);
  iDataP=length(iData[1,]);
  Data=NULL;
  for(i in 1:iDataNum)
  {
    for(j in 1:iDataP)
      Data[(i-1)*iDataP+j]=iData[i,j];
  }
  long_score <- vector("double",3*iDataP*(iDataP-1)/2)
  result <- .C("equSA_sub",
     as.numeric(Data),
     as.integer(iMaxNei),
     as.integer(iDataNum),
     as.integer(iDataP),
     as.numeric(ALPHA1),
     as.integer(GRID),
     as.integer(GRID),
     as.integer(iteration),
     as.numeric(long_score)
  )
  score <-matrix(result[[9]], ncol=3, byrow=TRUE)
  return(score)
}



ContTran <- function(iData,total_iteration=5000,stepsize=0.05)
{
  iDataNum=length(iData[,1]);
  iDataP=length(iData[1,]);
  iData=data.frame(t(iData))
  Data=NULL;
  for(i in 1:iDataP)
  {
    for(j in 1:iDataNum)
      Data[(i-1)*iDataNum+j]=iData[i,j];
  }
  cont_ave <- vector("double",iDataP*iDataNum)
  result <- .C("allR",
     as.integer(Data),
     as.integer(iDataNum),
     as.integer(iDataP),
     as.integer(total_iteration),
     as.numeric(stepsize),
     as.numeric(cont_ave)
  )
  
  continuz <- matrix(result[[6]],ncol=iDataP,byrow=TRUE)
  return(continuz)
}

 pcorselR <- function(score,ALPHA2=0.05,GRID=2,iteration=100)
 {
   irow=score[,1];
   icol=score[,2];
   idatax=score[,3];
   leng=length(irow);
   aaa <- vector("double",5*leng)
   result <- .C("pcorsel1",
      as.integer(irow),
      as.integer(icol),
      as.numeric(idatax),
      as.integer(leng),
      as.integer(GRID),
      as.integer(GRID),
      as.integer(iteration),
      as.numeric(aaa)
   )
   
   aaa<-matrix(result[[8]], ncol=5, byrow=TRUE)
   qqq<-ALPHA2
   FFF<-aaa[aaa[,5]<qqq,]
   if(sum(abs(FFF))==0){print("ALPHA2 should be larger!"); qqqscore=max(aaa[,3]);
   }else if(length(FFF)==5){
     print("ALPHA2 should be larger! Only one pair is selected!");
     qqqscore <- FFF[3];
   }else{
     qqqscore<-FFF[1,3]
   }
   return(qqqscore)
 }
 
 
 Mulpval <- function(pvalue,ALPHA2=0.05,GRID=2,iteration=100)
 {
   z <- -qnorm(pvalue,lower.tail=TRUE)
   leng <- length(z)
   irow=1:leng;
   icol=irow;
   idatax=z;
   aaa <- vector("double",5*leng)
   result <- .C("pcorsel1",
                as.integer(irow),
                as.integer(icol),
                as.numeric(idatax),
                as.integer(leng),
                as.integer(GRID),
                as.integer(GRID),
                as.integer(iteration),
                as.numeric(aaa)
   )
   
   aaa<-matrix(result[[8]], ncol=5, byrow=TRUE)
   qqq<-ALPHA2
   FFF<-aaa[aaa[,5]<qqq,]
   if(sum(abs(FFF))==0){print("ALPHA2 should be larger!"); qqqscore=max(aaa[,3]);
   }else if(length(FFF)==5){
     print("ALPHA2 should be larger! Only one pair is selected!");
     qqqscore <- FFF[3];
   }else{
     qqqscore<-FFF[1,3]
   }
   return(pnorm(-qqqscore))
 }
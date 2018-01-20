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
 
   .C("equSA1",
      as.numeric(Data),
      as.integer(iMaxNei),
      as.integer(iDataNum),
      as.integer(iDataP),
      as.numeric(ALPHA1),
      as.numeric(ALPHA2),
      as.integer(GRID),
      as.integer(GRID),
      as.integer(iteration)
      )
   A <- matrix(scan("sim0.mat"),ncol=iDataP,byrow=TRUE)
   score <-matrix(scan("sim0.pcor.est.ini"), ncol=3, byrow=TRUE)
   result <- list()
   result$Adj = A
   result$score = score
   unlink("sim0.mat")
   unlink("sim0.pcor.est.ini")
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
  
  .C("equSA_sub",
     as.numeric(Data),
     as.integer(iMaxNei),
     as.integer(iDataNum),
     as.integer(iDataP),
     as.numeric(ALPHA1),
     as.integer(GRID),
     as.integer(GRID),
     as.integer(iteration)
  )
  score <-matrix(scan("sim0.pcor.est.ini"), ncol=3, byrow=TRUE)
  unlink("sim0.pcor.est.ini")
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
  
  .C("allR",
     as.integer(Data),
     as.integer(iDataNum),
     as.integer(iDataP),
     as.integer(total_iteration),
     as.numeric(stepsize)
  )
  
  continuz <- matrix(scan("continuzed.dat"),ncol=iDataP,byrow=TRUE)
  unlink("continuzed.dat")
  return(continuz)
}

 pcorselR <- function(score,ALPHA2=0.05,GRID=2,iteration=100)
 {
   irow=score[,1];
   icol=score[,2];
   idatax=score[,3];
   length=length(irow);
   .C("pcorsel1",
      as.integer(irow),
      as.integer(icol),
      as.numeric(idatax),
      as.integer(length),
      as.integer(GRID),
      as.integer(GRID),
      as.integer(iteration)
   )
   
   aaa<-matrix(scan("aaa.fdr"), ncol=5, byrow=TRUE)
   unlink("aaa.fdr")
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
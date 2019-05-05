

MNR <- function(x,y,family="gaussian",penalty='lasso',tune='bic',alpha1=0.1,alpha2=0.05,level=0.95){
  p <- dim(x)[2]
  ## MNR method ##
  fit_SIS <- SIS(x,y,family = family,penalty = penalty, tune=tune)
  penal_sel <- fit_SIS$ix
  
  Res <- equSAR(x,ALPHA1=alpha1,ALPHA2=alpha2)
  A1<-Res$Adj
  
  CI <- NULL;
  beta <- NULL;
  pvalue <- NULL;
  for(k in 1:p){
    neibor <- which(A1[k,]==1)
    x_sel_name <- unique(c(neibor,penal_sel,k))
    x_sel <- x[,match(x_sel_name,1:p)]
    colnames(x_sel) <- x_sel_name
    if(family=="gaussian"){
      fit_lm2 <- lm(y~x_sel)
      CI <- rbind(CI,confint(fit_lm2,level=level)[which(x_sel_name==k)+1,])
      beta <- c(beta,fit_lm2$coefficients[which(x_sel_name==k)+1])
      pvalue <- c(pvalue,summary(fit_lm2)$coefficients[which(x_sel_name ==k)+1,4])
    }else if(family=="binomial"){
      fit_lm2 <- glm(y~x_sel,family=binomial(link="logit"))
      CI <- rbind(CI,confint(fit_lm2,level=level)[which(x_sel_name==k)+1,])
      beta <- c(beta,fit_lm2$coefficients[which(x_sel_name==k)+1])
      pvalue <- c(pvalue,summary(fit_lm2)$coefficients[which(x_sel_name ==k)+1,4])
    }else if(family=="cox"){
      fit_lm2 <- coxph(y ~ x_sel)
      CI <- rbind(CI,confint(fit_lm2,level=level)[which(x_sel_name==k),])
      beta <- c(beta,fit_lm2$coefficients[which(x_sel_name==k)])
      pvalue <- c(pvalue,summary(fit_lm2)$coefficients[which(x_sel_name ==k),5])
    }else{
      stop("Only 'gaussian', 'binomial' and 'cox' models are provided.\n")
    }
    

  }
  return(list(CI=CI,coef=beta,pvalue=pvalue))
  # CI  ## Confidence Intervals for all coefficients
  # beta  ## estimated regression coefficients 
  # pvalue  ## pvalues
  
  
}



plearn.moral<- function(data, gaussian.index = NULL, binary.index = NULL, 
                        alpha1 = 0.1, alpha2 = 0.02, restrict = FALSE, score.only=FALSE) {
  if (is.null(data))
    stop("Data is missing!")
  DataP <- dim(data)[2]
  DataNum <- dim(data)[1]
  if (is.null(gaussian.index)||is.null(binary.index)){
    binary.index <- (1:DataP)[apply(data,2,function(s) all(s%%1==0))]
    gaussian.index <- setdiff(1:DataP,binary.index)
  }
  
  
  if(restrict){
    p1 <- pairwise_p_fast.restric(data, gaussian.index = gaussian.index, binary.index = binary.index)
  }else{
    p1 <- pairwise_p_fast(data, gaussian.index = gaussian.index, binary.index = binary.index)
  }
  
  ## subsample
  p1<-p1[order(p1[,3]), 1:3]
  N<-dim(p1)[1]
  ratio<-ceiling(N/100000)
  m<-floor(N/ratio)
  m0<-N-m*ratio
  s<-sample.int(ratio,m,replace=TRUE)
  for(i in 1:length(s)) s[i]<-s[i]+(i-1)*ratio
  if(m0>0){
    s0<-sample.int(m0,1)+length(s)*ratio
    s<-c(s,s0)
  }
  p1s<-p1[s,]
  qscore1 = pcorselR(score=p1s, ALPHA2=alpha1)
  
  #### sep: separator
  #### margmatrix: marginal graph
  sep <- matrix(0, DataP , DataP + 1)
  SS <- matrix(0, DataP, DataP)
  margmatrix <- matrix(0, DataP, DataP)
  
  
  for(i in 1:dim(p1)[1]) {
    k1 = p1[i,1]
    k2 = p1[i,2]
    un = p1[i,3]
    if (un > qscore1) {
      sep[k1,DataP+1] = sep[k1,DataP+1]+1
      sep[k1,sep[k1,DataP+1]]=k2
      #sep[k2,DataP+1]= sep[k2,DataP+1]+1
      #sep[k2,sep[k2,DataP+1]]=k1
      
      SS[k1,sep[k1,DataP+1]] = un
      #SS[k2,sep[k2,DataP+1]] = un
      
      
    }
  }
  
  #########################################reduce neigbor
  DataNum = dim(data)[1]
  MaxNei = DataNum/(log(DataNum))
  for(i in 1:DataP) {
    if(sep[i,DataP+1] > MaxNei) {
      m = sep[i, DataP+1]
      ss_order = sort(SS[i, 1:m])
      indx <- match(ss_order, SS[i, 1:m])
      S = sep[i, 1:m]
      km = m - MaxNei
      sep[i, 1:MaxNei] = S[indx[(km+1):(km+MaxNei)]]
      sep[i, DataP+1] = MaxNei
      
    }
  }
  
  ###########        marginal adjacency matrix  ########################
  for(i in 1:DataP){
    for(m in 1:sep[i, DataP+1]) {
      margmatrix[i, sep[i, m]] = 1
      margmatrix[sep[i, m], i] = 1
    }
  }
  
  
  ####################################################################
  ###################    detect v-structure  #########################
  #### if i-k-j but i and j are not connected, keep i-j
  
  
  
  moralmatrix <- v_structure_test(data, sep, margmatrix,
                                  gaussian.index = gaussian.index, binary.index = binary.index)
  
  # moralmatrix <- margmatrix
  ## keep married parents in sep2 (parents are connected)
  sep2 <- matrix(0, DataP, DataP+1)
  for(i in 1:DataP) {
    for(j in 1:DataP) {
      if (moralmatrix[i, j] == 1) {
        sep2[i, DataP+1] = sep2[i, DataP+1] + 1
        sep2[i, sep2[i, DataP+1]] = j
      }
    }
  }
  
  
  if(restrict){
    p2 <- contional_p_fast.restric(data, sep2, gaussian.index = gaussian.index, binary.index = binary.index)
  }else{
    p2 <- contional_p_fast(data, sep2, gaussian.index = gaussian.index, binary.index = binary.index)
  }
  score1 <- p2
  if(score.only){
    return(score1)
  }else{
    
    p2<-p2[order(p2[,3]), 1:3]
    N<-dim(p2)[1]
    ratio<-ceiling(N/100000)
    m<-floor(N/ratio)
    m0<-N-m*ratio
    s<-sample.int(ratio,m,replace=TRUE)
    for(i in 1:length(s)) s[i]<-s[i]+(i-1)*ratio
    if(m0>0){
      s0<-sample.int(m0,1)+length(s)*ratio
      s<-c(s,s0)
    }
    p2s<-p2[s,]
    
    qscore2 = pcorselR(score=p2s, ALPHA2=alpha2)
    
    
    E <- matrix(0, DataP, DataP)
    
    for(i in 1:dim(p2)[1]) {
      if(p2[i,3] > qscore2) {
        E[p2[i,1], p2[i,2]] = 1
        E[p2[i,2], p2[i,1]] = 1
      }
      
    }
    
    
    return (list(moral.matrix = E,score=score1))
  }
  
}





pairwise_p_fast <- function(xy, gaussian.index = NULL, poisson.index = NULL, binary.index = NULL) {
  DataP = dim(xy)[2]
  score_matrix = NULL
  for(i in 1:DataP) {
    for(j in 1:DataP) {
      if(i != j) {
        if (i %in% gaussian.index) {
          family = gaussian()
        } else if (i %in% poisson.index) {
          family = poisson()
        } else family = binomial()
        glm1 <- speedglm(xy[,i] ~ xy[,j], family = family )
        qscore <- (-1)*qnorm( as.numeric(as.character(summary(glm1)$coeff[2,4])))
        qscore[qscore == Inf] = 35.0 + rnorm(1, 0, 0.1)
        qscore[qscore == -Inf] = -35.0 + rnorm(1, 0, 0.1)
        score_matrix <- rbind(score_matrix, c(i, j, round(qscore, 300)))
      }
    }
  }
  return (pscore_1 = score_matrix)
}











####conditional screening

contional_p_fast <- function(xy, sep2,
                             gaussian.index = NULL, poisson.index = NULL, binary.index = NULL){
  partial_score <- NULL
  DataP = dim(xy)[2]
  xy <- as.data.frame(xy)
  colnames(xy) = paste("x", 1:DataP, sep = "")
  for(i in 1:DataP) {
    if(i %in% gaussian.index) {
      family = gaussian()
    } else if (i %in% poisson.index) {
      family = poisson()
    } else family = binomial()
    for(j in 1:DataP) {
      
      if(i!=j) {
        if (sep2[i,DataP+1] < sep2[j, DataP+1]){
          m = sep2[i, DataP+1]
          if(m > 0) {
            neig = sep2[i, 1:m]
          } else neig  = NULL
          
        } else {
          m = sep2[j, DataP + 1]
          if(m > 0) {
            neig = sep2[j, 1:m]
          } else neig  = NULL
          
        }
        neig <- setdiff(neig, c(i, j))
        bneig  = unique(c(j, neig))
        sneig = setdiff(bneig, j)
        bform <- as.formula(paste(paste(paste("x", i, sep = ""), "~", sep=""),
                                  paste(paste("x",bneig,sep=""),collapse="+")))   
        l1 <- speedglm(bform, data = xy, family = family, model = TRUE, y = TRUE)
        if(length(sneig) == 0) {
          result <- drop1(l1, data = xy, test = "Chisq")$"Pr(>Chi)"[2]
        } else {
          sform <- as.formula(paste(paste("scope", "~", sep=""),
                                    paste("x", j, sep=""))) 
          result <- drop1(l1, sform, data = xy, test = "Chisq")$"Pr(>Chi)"[2]
        }
        
        result[which(is.na(result))] <- 1
        score <- round(-1*qnorm(result),300)
        score[score == Inf] <- 35 + rnorm(1, 0, 0.1)
        score[score == -Inf] <- -35 + rnorm(1, 0, 0.1)
        
        partial_score <- rbind(partial_score, c(i, j, round(score, 300)))
        
      }
    }
  }
  return (pscore_2 = partial_score)
  
  
}




######reduce the size of spouse
v_structure_test <- function(xy, sep, margmatrix,
                             gaussian.index = NULL, poisson.index = NULL, binary.index = NULL) {
  
  xy <- as.data.frame(xy)
  DataP = dim(xy)[2]
  colnames(xy) = paste("x", 1:DataP, sep = "")
  moralmatrix = margmatrix
  
  ##### detect v-structures
  for(i in 1:DataP) {
    for(j in 1:DataP) {
      if(i != j & margmatrix[i, j] == 0) {
        if(sep[i, DataP+1] > 0 & sep[j, DataP+1] >0 ) {
          kindx <- intersect(sep[i, 1:sep[i, DataP+1]], sep[j, 1:sep[j, DataP+1]])
          if(length(kindx) > 0) {
            ####test whether v-structure exists
            if (i %in% gaussian.index) {
              family = gaussian()
            } else if (i %in% poisson.index) {
              family = poisson()
            } else family = binomial()
            for(kk in kindx) {
              bform <- as.formula(paste(paste(paste("x", i, sep = ""), "~", sep=""),
                                        paste(paste("x",c(j, kk),sep=""),collapse="+"))) 
              sform <- as.formula(paste(paste("scope", "~", sep=""),
                                        paste("x", j, sep=""))) 
              g1 <- speedglm(bform, data = xy, family = family, model = TRUE, y = TRUE)
              result <- drop1(g1, sform, data = xy, test = "Chisq")$"Pr(>Chi)"[2]
              #cat("spouse: ", i, j, "child:", kk,result, "\n")
              if (result < 0.1) {
                #cat("spouse: ", i, j, "child:", kk,result, "\n")
                moralmatrix[i, j] = 1
              }
            }
            
            
          }
        }
      }
    }
  }
  return (moralmatrix)
}




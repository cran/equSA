

#######resolving Markov Blanket
Triconflict <- function(xy, E, i, j, Trij, family = c("gaussian", "binomial", "poisson")) {
  this.call = match.call()
  family = match.arg(family)
  if (is.null(xy))
    stop("Data is missing!")
  fit = switch(family, gaussian = triconflictglm(xy, E, i, j, Trij, "gaussian"), binomial = 
                 triconflictglm(xy, E, i, j, Trij, "binomial"), 
               poisson = triconflictglm(xy, E, i, j, Trij, "poisson"))
  fit$call = this.call
  #class(fit) = c(class(fit), "Triconflict")
  return(fit)  
}    

triconflictglm <- function(xy, E, i, j, Trij, family) {
  xy <- as.data.frame(xy)
  DataP <- dim(xy)[2]
  colnames(xy) <- paste("x", 1:DataP, sep = "")
  Bdi = which(E[i, ] == 1)
  Bdj = which(E[j, ] == 1)
  Bi = setdiff(Bdi, union(Trij, j))
  Bj = setdiff(Bdj, union(Trij, i))
  if (length(Bi) > length(Bj)) {
    B = Bj
  } else {B = Bi}
  #cat("B", B, "\n")
  ## list all strict subset
  subset = NULL
  for(s in 0:(length(Trij) - 1)) {
    subset = c(subset, combn(Trij, s, simplify = FALSE))}
  
  for (m in 1:length(subset)) {
    subsetm = subset[[m]]
    z = union(B, subsetm)
    #cat("z", z, "\n")
    if (length(z) > 1) {
      reg1 <- paste("x", z, sep = "")
      reg2 <-  paste("x", c(j, z), sep = "")
      
      #cat("reg1", reg1, "\n")
      #cat("reg2", reg2, "\n")
      response <- paste("x", i, sep ="")
      response <- paste(response, "~", sep="")
      fom1 <- as.formula(paste(response, paste(paste(reg1, collapse = "+"), sep = "+")))
      fom2 <- as.formula(paste(response, paste(paste(reg2, collapse = "+"),sep = "+")))
      
      l1 <- glm(fom1, data = xy, family )
      l2 <- glm(fom2, data = xy, family )
      
      rij = anova(l1, l2, test = "Chisq")$"Pr(>Chi)"[2] } else {
        
        
        glm1 <- glm(xy[,i] ~ xy[,j], data = xy, family )
        rij <- summary(glm1)$coeff[2,4]
        
      }
    #cat("rij", rij, "\n")
    return (list(i=i, j=j, z=z, rij = rij))
    
    ## reachable nodes by w
    w = setdiff(Trij, subsetm)
    g_ij = setdiff(1:DataP, c(i, j))
    
    
    if (length(w) >= 1) {
      w_reach = NULL
      for(wr in w) {
        dis = distances(graph_from_adjacency_matrix(E, mode = "undirected"),
                        v = wr, to = setdiff(g_ij, wr))
        w_reach = c(w_reach, setdiff(g_ij, wr)[which(dis < Inf)])
      }
      w_reach = unique(w_reach)
    } else w_reach = NULL
    D = intersect(B, w_reach)
    #cat("D", D, "\n")
    B_prime = setdiff(B, D)
    #cat("B_prime", B_prime, "\n")
    ### susbet of D
    Dsubset = NULL
    if (length(D) > 0) {
      for(Ds in 0:(length(D) - 1)) {
        Dsubset = c(Dsubset, combn(D, Ds, simplify = FALSE))}
      
      
      for (Dm in 1:length(Dsubset)) {
        Dsubsetm = Dsubset[[Dm]] 
        
        z = Reduce(union, list(B_prime, Dsubsetm, subsetm))
        if (length(z) > 1) {
          reg1 <- paste("x", z, sep = "")
          reg2 <-  paste("x", c(j, z), sep = "")
          
          
          response <- paste("x", i, sep ="")
          response <- paste(response, "~", sep="")
          fom1 <- as.formula(paste(response, paste(paste(reg1, collapse = "+"), sep = "+")))
          fom2 <- as.formula(paste(response, paste(paste(reg2, collapse = "+"),sep = "+")))
          
          l1 <- glm(fom1, data = xy, family )
          l2 <- glm(fom2, data = xy, family )
          
          
          rij = anova(l1, l2, test = "Chisq")$"Pr(>Chi)"[2] 
        } else {
          
          
          glm1 <- glm(xy[,i] ~ xy[,j], data = xy)
          rij <- summary(glm1)$coeff[2,4]
          
        }
        return(list(i=i, j=j, z=z, rij = rij))
        
        
        
        
      } 
    }
    
  }  
  
}  



############################Resolve Markov blanket with conflict sets
####input E: moral graph
####variable type
####output G: partially oriented DAG
####spouse index: spouse_index
ResolveMB <- function(E, xy, gaussian_index = NULL, poisson_index = NULL, binomial_index = NULL, alpha = 0.1) {
  G = E
  spouse_ind <- NULL
  Cset = NULL
  DataP = dim(E)[1]
  for(i in 1:DataP) {
    for(j in 1:DataP) {
      
      if (E[i, j] == 1) {
        Trij = intersect(which(E[i, ] == 1), which(E[j, ] == 1))
        
        if (length(Trij) > 0) {
          
          if(i %in% poisson_index) {
            sep = Triconflict(xy, E, i, j, Trij, family = "poisson")
          } else if (i %in% binomial_index) {
            sep = Triconflict(xy, E, i, j, Trij, family = "binomial")
          } else{
            sep = Triconflict(xy, E, i, j, Trij, family = "gaussian")
          }
          
          if(is.null(sep) == FALSE & sep$rij > alpha) {
            spouse_ind = rbind(spouse_ind, c(i, j))
            z = setdiff(Trij, sep$z)
            for (zm in z) {
              Cset = rbind(Cset, c(i, j, zm))
              
            }
          }    
        }
      }
      
    }
  }
  
  G[spouse_ind] = 0
  if(length(Cset) > 0) {
    for(i in 1:dim(Cset)[1]) {
      if (G[Cset[i, 1], Cset[i, 3]] == 1 & G[Cset[i, 2], Cset[i, 3]] == 1 &
            G[Cset[i, 3], Cset[i, 1]] == 1 & G[Cset[i, 3], Cset[i, 2]] == 1) {
        #cat("i, j, =0", Cset[i, 3], Cset[i, 1], "\n")
        #cat("i, j, =0", Cset[i, 3], Cset[i, 2], "\n")
        G[Cset[i, 3], Cset[i, 1]] = 0
        G[Cset[i, 3], Cset[i, 2]] = 0
      }
    }
  }
  
  return (list(PDAG = G, spouse_index = spouse_ind))
}


######pairwise screening


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
              if (result < 0.3) {
                #cat("spouse: ", i, j, "child:", kk, "\n")
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



p_learning<- function(data, gaussian.index = NULL, binary.index = NULL, poisson.index = NULL, 
                        alpha1 = 0.1, alpha2 = 0.02, alpha3 = 0.02) {
  if (is.null(data))
    stop("Data is missing!")
  DataP <- dim(data)[2]
  DataNum <- dim(data)[1]
  
  p1 = pairwise_p_fast(data, gaussian.index = gaussian.index, binary.index = binary.index,
                       poisson.index = poisson.index)
  
  qscore1 = pcorselR(p1, alpha1)
  
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
  MaxNei = DataNum/(5*log(DataNum))
  for(i in 1:DataP) {
    if(sep[i,DataP+1] > MaxNei) {
      m = sep[i, DataP+1]
      ss_order = sort(SS[i, 1:m])
      indx <- match(ss_order, SS[i, 1:m])
      S = sep[i, 1:m]
      k = m - MaxNei
      sep[i, 1:MaxNei] = S[indx[(k+1):(k+MaxNei)]]
      sep[i, DataP+1] = MaxNei
      
    }
  }
  ######################################################################
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
                                  gaussian.index = gaussian.index, binary.index = binary.index,
                                  poisson.index = poisson.index)
  
  
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
  
  
  p2 = contional_p_fast(data, sep2, gaussian.index = gaussian.index, binary.index = binary.index,
                        poisson.index = poisson.index)
  
  qscore2 = pcorselR(p2, alpha2)
  
  
  
  
  
  E <- matrix(0, DataP, DataP)
  
  for(i in 1:dim(p2)[1]) {
    if(p2[i,3] > qscore2) {
      #cat("i, j", p2[i, 1], p2[i,2], p2[i, 3], "\n")
      E[p2[i,1], p2[i,2]] = 1
      E[p2[i,2], p2[i,1]] = 1
    }
    
  }
  
  
  
  
  G = ResolveMB(E, data, gaussian_index = gaussian.index, 
                poisson_index = poisson.index, binomial_index = binary.index, alpha = alpha3)$PDAG
  
  return (list(PDAG = G))
  
}

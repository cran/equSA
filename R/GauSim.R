
## Main function
GauSim = function(n, p, graph = "AR(2)", v = NULL, u = NULL, g = NULL, prob = NULL, vis = FALSE, verbose = TRUE){	
  gcinfo(FALSE)
  if(verbose) cat("Generating data from the Gaussian distribution with the", graph,"graph structure....\n")
  if(is.null(g)){
    g = 1
    if(graph == "hub" || graph == "cluster"){
      if(p > 40)	g = ceiling(p/20)
      if(p <= 40) g = 2
    }
  }
  
  if(graph == "random"){
    if(is.null(prob))	prob = min(1, 3/p)
    prob = sqrt(prob/2)*(prob<0.5)+(1-sqrt(0.5-0.5*prob))*(prob>=0.5)
  }
  
  if(graph == "cluster"){
    if(is.null(prob)){
      if(p/g > 30)	prob = 0.3
      if(p/g <= 30)	prob = min(1,6*g/p)
    }
    prob = sqrt(prob/2)*(prob<0.5)+(1-sqrt(0.5-0.5*prob))*(prob>=0.5)
  }  
  
  
  # parition variables into groups
  g.large = p%%g
  g.small = g - g.large
  n.small = floor(p/g)
  n.large = n.small+1
  g.list = c(rep(n.small,g.small),rep(n.large,g.large))
  g.ind = rep(c(1:g),g.list)
  rm(g.large,g.small,n.small,n.large,g.list)
  gc()
  
  # build the graph structure
  theta = matrix(0,p,p);
  C = matrix(0,p,p);
  if(graph == "AR(2)"){
    if(is.null(u)) u = 0.1
    if(is.null(v)) v = 0.3
    for(i in 1:p)
    {
      for(j in 1:p)
      {
        if(abs(j-i)==1) {C[i,j]=0.5; theta[i,j]=1}
        else if(abs(j-i)==2) {C[i,j]=0.25; theta[i,j]=1}
      }
    }
    diag(C) = 0
    omega = C
  }
  if(graph == "cluster"){
    if(is.null(u)) u = 0.1
    if(is.null(v)) v = 0.3
    for(i in 1:g){
      tmp = which(g.ind==i)
      tmp2 = matrix(runif(length(tmp)^2,0,0.5),length(tmp),length(tmp))
      tmp2 = tmp2 + t(tmp2)		 	
      theta[tmp,tmp][tmp2<prob] = 1
      diag(theta) = 0
      omega = theta*v
      rm(tmp,tmp2)
      gc()
    }
  }
  if(graph == "hub"){
    if(is.null(u)) u = 0.1
    if(is.null(v)) v = 0.3
    for(i in 1:g){
      tmp = which(g.ind==i)
      theta[tmp[1],tmp] = 1
      theta[tmp,tmp[1]] = 1
      diag(theta) = 0
      omega = theta*v
      rm(tmp)
      gc()
    }
  }
  if(graph == "random"){
    if(is.null(u)) u = 0.1
    if(is.null(v)) v = 0.3
    
    tmp = matrix(runif(p^2,0,0.5),p,p)
    tmp = tmp + t(tmp)
    theta[tmp < prob] = 1
    diag(theta) = 0
    omega = theta*v
    #theta[tmp >= tprob] = 0
    rm(tmp)
    gc()
  }
  
  if(graph == "scale-free"){
    if(is.null(u)) u = 0.1
    if(is.null(v)) v = 0.3
    g <- igraph::barabasi.game(n=p, power=0.01, zero.appeal=p, directed=F)
    theta <- as.matrix(igraph::get.adjacency(g, type="both"))
    diag(theta) = 0
    omega = theta*v
  }
  
  
  # make omega positive definite and standardized
  diag(omega) = abs(min(eigen(omega)$values)) + 0.1 + u
  sigma = solve(omega)
  mu <- rep(0,p)
  data <- rmvnorm(n,mu,sigma)
  sim = list(data = data, sigma = sigma, theta = theta)
  return(sim)
}



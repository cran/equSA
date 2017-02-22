
ContSim <- function(n, p,  v = NULL, u = NULL, g = NULL, prob = NULL, vis = FALSE, verbose = TRUE, graph.type="AR(2)",
                    k=3.30, lambda=515, omega=0.003,lower.tail = TRUE, log.p = FALSE)
{

B <- SimGraph( p=p, graph = graph.type, v = v, u = u, g = g, prob = prob, vis = vis, verbose = verbose)
mu <- rep(0,p)
sigma <- B$sigma
x <- rmvnorm(n,mu,sigma)


# generate count data y following zero-inflated negative-binomial distribution 
y <- x
for(j in 1:p)
{
  dat <- x[,j]
  mu_v <- mean(dat)
  sd_v <- sd(dat)
  p_v <- pnorm(dat,mu_v,sd_v)
  y[,j] <- qzinb(p = p_v, k=k, lambda=lambda, omega=omega, lower.tail = lower.tail, log.p = log.p)
}
result <- list()
result$data = y
result$Adj = B$theta

return(result)
}

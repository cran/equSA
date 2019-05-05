
SimGraDat <- function(n = 200, p = 100, type = "band", rate = 0.1) {
  if (type != "band") {
    stop("Only 'band' type is provided !")
  }
  C <- diag(1, p)
  for (i in 1:p) {
    for (j in 1:p) {
      if (abs(j - i) == 1) 
        C[i, j] = 0.5 else if (abs(j - i) == 2) 
          C[i, j] = 0.25
    }
  }
  A <- diag(0, p)
  sigma <- solve(C)
  mu <- rep(0, p)
  sigma2 <- sigma
  for (i in 1:p) {
    for (j in 1:p) {
      if (sigma[i, j] <= sigma[j, i]) {
        sigma2[i, j] <- sigma[j, i]
      } else {
        sigma2[i, j] <- sigma[i, j]
      }
      if (C[i, j] != 0 && i != j) 
        A[i, j] <- 1
    }
  }
  data <- rmvnorm(n, mu, sigma2)
  data_true <- data
  length(data)
  sec <- sample(1:length(data), size = rate * length(data))
  data[sec] <- NA
  result <- list(data = data, A = A)
  return(result)
}

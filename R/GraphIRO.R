GraphIRO <- function(data, A, alpha1 = 0.05, alpha2 = 0.05, alpha3 = 0.05, iteration = 30, warm = 10)
{
  p <- dim(data)[2]
  leng <- iteration - warm
  miss <- which(is.na(data), arr.ind = TRUE)
  data_result <- data
  for(i in 1:length(miss[,1]))
  {
    data_result[miss[i,1],miss[i,2]] <- median(na.omit(data[,miss[i,2]]))
  }
  U <- NULL
  for (j in 1:iteration)
  {
    GraRes <- equSAR(data_result, ALPHA1 = alpha1, ALPHA2 = alpha2)
    U1 <- GraRes$score
    A2 <- GraRes$Adj
    U <- cbind(U, U1[, 3])
    for (i in 1:length(miss[, 1]))
    {
      dep <- which(A2[miss[i, 2], ] == 1)
      if (length(dep) > 0)
      {
        combine <- data.frame(data_result[, miss[i, 2]], data_result[, dep])
        cov <- cov(combine)
        mu1 <- mean(data_result[, miss[i, 2]])
        mu2 <- apply(data_result[, dep, drop = FALSE], 2, mean)
        mu <- mu1 + cov[-1, 1] %*% solve(cov[-1, -1]) %*% t(data_result[miss[i, 1], dep, drop = FALSE] - mu2)
        sigma <- cov[1, 1] - cov[-1, 1] %*% solve(cov[-1, -1]) %*% cov[1, -1]
        data_result[miss[i, 1], miss[i, 2]] <- rnorm(1, mu, sqrt(sigma))
      }
      else
      {
        next
      }
    }
  }
  N <- p * (p - 1)/2
  ratio <- ceiling(N/1e+05)
  U <- cbind(U1[, 1:2], U)
  z <- -apply(U[, -(1:(iteration - leng + 2)),drop=FALSE], 1, sum)/leng
  q <- pnorm(-abs(z), log.p = TRUE)
  q <- q + log(2)
  s <- qnorm(q, log.p = TRUE)
  s <- (-1) * s
  U <- cbind(U[, 1:2], s)
  fit <- summary(U[, 3])
  cut1 <- seq(fit[1], fit[2], length = 10)
  cut2 <- seq(fit[2], fit[3], length = 10)
  cut3 <- seq(fit[3], fit[5], length = 10)
  cut4 <- seq(fit[5], fit[6], length = 30)
  cut <- c(cut1, cut2, cut3, cut4)
  precision <- recall <- NULL
  w3 <- NULL
  for (i in 1:length(cut))
  {
    A2 <- diag(0, p)
    aa <- U[cut[i] - U[, 3] < 0.001, ,drop = FALSE]
    for (j in 1:nrow(aa))
    {
      A2[aa[j, 1], aa[j, 2]] <- 1
      A2[aa[j, 2], aa[j, 1]] <- 1
    }
    B <- A + A2
    L <- length(which(B == 2))
    N <- length(which(B == 0))
    precision <- 1 * L/sum(A2)
    recall <- 1 * L/sum(A)
    w3 <- rbind(w3, c(recall, precision))
  }
  RecPre <- w3
  colnames(RecPre) <- c("Recall", "Precision")
  psiscore <- U
  U <- U[order(U[, 3]), 1:3]
  N <- p * (p - 1)/2
  ratio <- ceiling(N/1e+05)
  m <- floor(N/ratio)
  m0 <- N - m * ratio
  s <- sample.int(ratio, m, replace = TRUE)
  for (i in 1:length(s)) s[i] <- s[i] + (i - 1) * ratio
  if (m0 > 0)
  {
    s0 <- sample.int(m0, 1) + length(s) * ratio
    s <- c(s, s0)
  }
  Us <- U[s, ]
  qqqscore <- pcorselR(Us, ALPHA2 = alpha3)
  U <- psiscore
  U <- U[U[, 3] >= qqqscore, ]
  Adj <- matrix(rep(0, p * p), ncol = p)
  for (i in 1:nrow(U))
  {
    k1 <- U[i, 1]
    k2 <- U[i, 2]
    Adj[k1, k2] <- 1
    Adj[k2, k1] <- 1
  }
  result <- list()
  result$Adj <- Adj
  result$RecPre <- RecPre
  return(result)
}

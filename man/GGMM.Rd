\name{GGMM}
\alias{GGMM}
\title{ Learning high-dimensional Gaussian Graphical Models with Heterogeneous Data.}
\description{Gaussian Graphical Mixture Models for learning a single high-dimensional network structure from heterogeneous dataset.}

\usage{
GGMM(data, A, M, alpha1 = 0.1, alpha2 = 0.05, alpha3 = 0.05, iteration = 30, warm = 20)
}

\arguments{
  \item{ data }{\eqn{n}x\eqn{p} mixture Gaussian distributed dataset.}
  \item{ A }{\eqn{p}x\eqn{p} true adjacency matrix for evaluating the performance.}
  \item{ M }{ The number of heterogeneous groups.}
\item{ alpha1  }{ The significance level of correlation screening in the \eqn{\psi}-learning algorithm, see R package \pkg{equSA} for detail. In general, a high significance level of correlation screening will lead to 
 a slightly large separator set, which reduces the risk of missing important variables in 
 the conditioning set. In general, including a few false variables in the conditioning set will not hurt much the 
 accuracy of the \eqn{\psi}-partial correlation coefficient, the default value is 0.1.}
  \item{ alpha2  }{ The significance level of \eqn{\psi}-partial correlation coefficient screening for estimating the adjacency matrix, see \pkg{equSA}, the default value is 0.05.}
\item{ alpha3  }{ The significance level of integrative \eqn{\psi}-partial correlation coefficient screening for estimating the adjacency matrix of GGMM method, the default value is 0.05.}
\item{ iteration  }{ The number of total iterations, the default value is 30.}
\item{ warm }{ The number of burn-in iterations, the default value is 20.}
}

\value{
   \item{ RecPre }{ The output of Recall and Precision values of our proposed method.}
  \item{ Adj }{  \eqn{p}x\eqn{p} Estimated adjacency matrix.}
  \item{ label }{ The estimated group indices for each observation.}
  \item{ BIC }{ The BIC scores for determining the number of groups \eqn{M}.}
}


\author{ Bochao Jia\email{jbc409@gmail.com} and Faming Liang}

\examples{
 \dontshow{
set.seed(10)
library(equSA)
result <- SimHetDat(n = 50, p = 10, M = 2, mu = 0.2, type = "band")
Est <- GGMM(result$data, result$A, alpha1 = 0.5, alpha2 = 0.5, alpha3 = 0.5, M = 2, iteration = 1, warm = 0)
}
 \donttest{
library(equSA)
result <- SimHetDat(n = 100, p = 200, M = 3, mu = 0.5, type = "band")
Est <- GGMM(result$data, result$A, M = 3, iteration = 30, warm = 20)
## plot network by our estimated adjacency matrix.
plotGraph(Est$Adj)
## plot the Recall-Precision curve
plot(Est$RecPre[,1], Est$RecPre[,2], type="l", xlab="Recall", ylab="Precision")  
}

}

\references{
Liang, F., Song, Q. and Qiu, P. (2015). An Equivalent Measure of Partial Correlation Coefficients for High Dimensional Gaussian Graphical Models. J. Amer. Statist. Assoc., 110, 1248-1265.

Liang, F. and Zhang, J. (2008) Estimating FDR under general dependence using stochastic approximation.  Biometrika,  95(4), 961-977.

Liang, F., Jia, B., Xue, J., Li, Q., and Luo, Y. (2018). An Imputation Regularized Optimization Algorithm for High-Dimensional Missing Data Problems and Beyond. Submitted to Journal of the Royal Statistical Society Series B.

Jia, B. and Liang, F. (2018). Learning Gene Regulatory Networks with High-Dimensional Heterogeneous Data. Accept by ICSA Springer Book. 


}

\keyword{GGMM}

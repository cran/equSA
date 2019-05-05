\name{MNR}
\alias{MNR}
\title{ Markov Neighborhood Regression for High-Dimensional Inference.}
\description{Construct confidence intervals and assess p-values in high-dimensional linear and generalized linear models.}

\usage{
MNR(x,y,family='gaussian',penalty='lasso',tune='bic',alpha1=0.1,alpha2=0.05,level=0.95)
}

\arguments{
  \item{ x }{The design matrix, of dimensions \eqn{n}x\eqn{p}, without an intercept. Each row is an observation vector.}
  \item{ y }{ The response vector of dimension \eqn{n}x\eqn{1}. Quantitative for family='gaussian', binary (0-1) for family='binomial'. For family='cox', y should be an object of class \code{Surv}, as provided by the function \code{Surv()} in the package \pkg{survival}.}
  \item{ family  }{Response type (see above).}
  \item{ penalty }{
  The penalty to be applied in the regularized likelihood subproblems. 'lasso' (the default), 'MCP', or 'SCAD' are provided. See package \pkg{SIS} for detail.}
  \item{ tune }{
  Method for tuning the regularization parameter of the penalized likelihood subproblems and of the final model selected by (I)SIS. Options include tune='bic', tune='ebic', tune='aic', and tune='cv'.}
\item{ alpha1  }{ The significance level of correlation screening in the \eqn{\psi}-learning algorithm, see R package \pkg{equSA} for detail. In general, a high significance level of correlation screening will lead to 
 a slightly large separator set, which reduces the risk of missing important variables in 
 the conditioning set. In general, including a few false variables in the conditioning set will not hurt much the 
 accuracy of the \eqn{\psi}-partial correlation coefficient, the default value is 0.1.}
  \item{ alpha2  }{ The significance level of \eqn{\psi}-partial correlation coefficient screening for estimating the adjacency matrix, see \pkg{equSA}, the default value is 0.05.}
   \item{ level }{the confidence level required, the default value is 0.95}
  
}




\value{
   \item{ CI }{ Estimated confidence intervals for all coefficients.}
  \item{ coef }{  \eqn{p}x\eqn{1} estimated regression coefficients for all variables.}
    \item{ pvalue }{  \eqn{p}x\eqn{1} estimated p-values for all variables.}
%% ...
}

\author{ Bochao Jia\email{jbc409@gmail.com} and Faming Liang}

\examples{
\dontshow{
library(equSA)
p <- 30
coef_true <- rep(0,p)
coef_true[1:5] <- c(2,4,-3,-5,10)
coef <- c(1,coef_true)
data <- SimMNR(n = 100, p = 30, coef = coef, family = "gaussian")
MNR(data$x, data$y, family = "gaussian")

}
 \donttest{
library(equSA)
p <- 500
coef_true <- rep(0,p)
coef_true[1:5] <- c(2,4,-3,-5,10)
coef <- c(1,coef_true)
data <- SimMNR(n = 200, p = 500, coef = coef, family = "gaussian")
MNR(data$x, data$y, family = "gaussian")

}

}

\references{

Liang, F., Xue, J. and Jia, B. (2018). Markov Neighborhood Regression for High-Dimensional Inference. Submitted to J. Amer. Statist. Assoc.

}

\keyword{MNR}

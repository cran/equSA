\name{SimMNR}
\alias{SimMNR}
\title{ Simulate Data for high-dimensional inference}
\description{Simulate data with graphical structure for generalized regression, which can be used in \code{MNR(x,y,...)} for constructing confidence intervals and assessing p-values.}

\usage{
SimMNR(n, p, coef, family="gaussian")
}

\arguments{
  \item{ n }{ Number of observations.}
  \item{ p }{ Number of variables.}
  \item{ coef }{A \eqn{p+1}x\eqn{1} vector. The first value denotes the intercept term and other \eqn{p} values denote the true regression coefficients for \eqn{p} variables.}
  \item{ family }{ Quantitative for family='gaussian' (default), binary (0-1) for family='binomial'. Survival data for family='cox'.}
}

\details{
We generate \eqn{p} variables from the following precision matrix, which is often been called "band" structure or "AR(2)" structure.
\deqn{
      C_{i,j}=\left\{\begin{array}{ll}
                     0.5,&\textrm{if $\left| j-i \right|=1, i=2,...,(p-1),$}\\
                     0.25,&\textrm{if $\left| j-i \right|=2, i=3,...,(p-2),$}\\
                     1,&\textrm{if $i=j, i=1,...,p,$}\\
                    0,&\textrm{otherwise.}
                \end{array}\right.
                }
}

\value{
   \item{ x }{ Simulated data in a \emph{n}x\emph{p} design matrix, without an intercept.}
  \item{ y }{  The response vector of dimension \eqn{n}x\eqn{1}. Quantitative for family='gaussian', binary (0-1) for family='binomial'. For family='cox', y should be an object of class \code{Surv}, as provided by the function \code{Surv()} in the package \pkg{survival}.}
  \item{A}{The true adjacency matrix of variables in the design matrix \eqn{x}.}
%% ...
}

\author{ Bochao Jia\email{jbc409@gmail.com} and Faming Liang}

\examples{
library(equSA)
p <- 200
coef_true <- rep(0,p)
coef_true[1:5] <- runif(5,3,5)
coef <- c(1,coef_true)
data <- SimMNR(n = 100, p = 200, coef = coef, family = "cox")
}

\references{

Liang, F., Xue, J. and Jia, B. (2018). Markov Neighborhood Regression for High-Dimensional Inference. Submitted to J. Amer. Statist. Assoc.

}

\keyword{SimMNR}

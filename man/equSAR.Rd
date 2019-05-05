\name{equSAR}
\alias{equSAR}
\title{ An equvalent mearsure of partial correlation coeffients}
\description{Infer networks from Gaussian data by \eqn{\psi}-learning algorithm.}

\usage{
equSAR(iData,iMaxNei,ALPHA1=0.05,ALPHA2=0.05,GRID=2,iteration=100)
}

\arguments{
  \item{ iData }{ a \eqn{n}x\eqn{p} data matrix.}
  \item{ iMaxNei  }{ Neiborhood size in correlation screening step, default to \eqn{n/log(n)}, where \eqn{n} is the number of observation.}
  \item{ ALPHA1  }{ The significance level of correlation screening. In general, a high significance level of correlation screening will lead to 
 a slightly large separator set \eqn{S_{ij}}, which reduces the risk of missing some important variables in 
 the conditioning set. Including a few false variables in the conditioning set will not hurt much the 
 accuracy of the \eqn{\psi}-partial correlation coefficient.}
  \item{ ALPHA2  }{ The significance level of \eqn{\psi} screening.}
  \item{GRID}{The number of components for the \eqn{\psi} scores. The default value is 2.}
    \item{iteration}{Number of iterations for screening. The default value is 100.}
}
\details{
This is the main function of the package that fit the Gaussian Graphical Models and obtain the \eqn{\psi} scores and adjacency matrix.
}
\value{
A list of two elements:
  \item{Adj}{\eqn{p}x\eqn{p} adjacency matrix of the generated graph.}
    \item{score}{Estimated \eqn{\psi} score matrix which has 3 columns. The first two columns denote the pair indices of variables \eqn{i} and \eqn{j} and the last column denote the calculated \eqn{\psi} scores for this pair.}
%% ...
}

\author{ Bochao Jia and Faming Liang}

\examples{
\dontshow{
library(equSA)
data <- GauSim(100,20)$data
equSAR(data,GRID=2,iteration=20)
}
\donttest{
library(equSA)
data <- GauSim(100,100)$data
equSAR(data)
}

}

\references{
Liang, F., Song, Q. and Qiu, P. (2015). An Equivalent Measure of Partial Correlation Coefficients for High Dimensional Gaussian Graphical Models. J. Amer. Statist. Assoc., 110, 1248-1265.

Liang, F. and Zhang, J. (2008) Estimating FDR under general dependence using stochastic approximation.  Biometrika,  95(4), 961-977.
}

\keyword{equSAR}

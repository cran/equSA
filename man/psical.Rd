\name{psical}
\alias{psical}
\title{ An calculation of \eqn{\psi} scores.}
\description{To compute an equvalent mearsure of partial correlation coeffients called \eqn{\psi} scores.}

\usage{
psical(iData,iMaxNei=as.integer(iDataNum/log(iDataNum)),
ALPHA1=0.05,GRID=2,iteration=100)
}

\arguments{
  \item{ iData }{ a \emph{n}x\emph{p} data matrix.}
  \item{ iMaxNei  }{ Neiborhood size in correlation screening step, default to n/log(n).}
  \item{ ALPHA1  }{ The significance level of correlation screening. In general, a high significance level of correlation screening will lead to 
 a slightly large separator set \eqn{S_{ij}}, which reduces the risk of missing some important variables in 
 the conditioning set. Including a few false variables in the conditioning set will not hurt much the 
 accuracy of the \eqn{\psi}-partial correlation coefficient.}
   \item{GRID}{The number of components for the corrlation scores. The default value is 2.}
\item{iteration}{Number of iterations for screening. The default value is 100.}
}
\details{
This is the function to calculate \eqn{\psi} scores and can be used in combining or detecting difference of two networks.
}
\value{
  \item{score}{Estimated \eqn{\psi} score matrix which has 3 columns. The first two columns denote the pair indices of variables i and j and the last column denote the calculated \eqn{\psi} scores for this pair.}
%% ...
}

\author{ Bochao Jia, Faming liang\email{faliang@ufl.edu}}

\examples{

\donttest{
library(equSA)
data(SR0)
subset <- SR0
psical(subset)
}
         
}

\references{
Liang, F., Song, Q. and Qiu, P. (2015). An Equivalent Measure of Partial Correlation Coefficients for High Dimensional Gaussian Graphical Models. J. Amer. Statist. Assoc., 110, 1248-1265.

Liang, F. and Zhang, J. (2008) Estimating FDR under general dependence using stochastic approximation.  Biometrika,  95(4), 961-977.
}

\keyword{psical}

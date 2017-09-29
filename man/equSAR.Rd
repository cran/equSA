\name{equSAR}
\alias{equSAR}
\title{ An equvalent mearsure of partial correlation coeffients}
\description{Infer networks from Gaussian data using our proposed \eqn{\psi} -learning algorithm.}

\usage{
equSAR(iData,iMaxNei=as.integer(iDataNum/log(iDataNum)),
ALPHA1=0.05,ALPHA2=0.05)
}

\arguments{
  \item{ iData }{ a \emph{n}x\emph{p} data matrix.}
  \item{ iMaxNei  }{ Neiborhood size in correlation screening step, default to n/log(n).}
  \item{ ALPHA1  }{ The significance level of correlation screening. In general, a high significance level of correlation screening will lead to 
 a slightly large separator set S_{ij}, which reduces the risk of missing some important variables in 
 the conditioning set. Including a few false variables in the conditioning set will not hurt much the 
 accuracy of the \eqn{\psi}-partial correlation coefficient.}
  \item{ ALPHA2  }{ The significance level of \eqn{\psi} screening.}
}
\details{
This is the main function of the package that fit the Gaussian Graphical Models and obtain the \eqn{\psi} scores and adjacency matrix.
}
\value{
A list of two elements:
  \item{Adj}{\emph{p}x\emph{p} adjacency matrix of the generated graph.}
    \item{score}{Estimated \eqn{\psi} score matrix which has 3 columns. The first two columns denote the pair indices of variables i and j and the last column denote the calculated \eqn{\psi} scores for this pair.}
%% ...
}

\author{ Bochao Jia and Faming Liang\email{faliang@ufl.edu}}

\examples{
library(equSA)
data(SR0)
subset <- SR0[1:100,1:10]
equSAR(subset)
      \dontshow{
            library(equSA)
            data(SR0)
            subset <- SR0[1:100,1:10]
            equSAR(subset)
      }
      \donttest{
            #library(equSA)
            #data(SR0)
            #equSAR(SR0)
      }
}

\references{
Liang, F., Song, Q. and Qiu, P. (2015). An Equivalent Measure of Partial Correlation Coefficients for High Dimensional Gaussian Graphical Models. J. Amer. Statist. Assoc., 110, 1248-1265.

Liang, F. and Zhang, J. (2008) Estimating FDR under general dependence using stochastic approximation.  Biometrika,  95(4), 961-977.
}

\keyword{equSAR}

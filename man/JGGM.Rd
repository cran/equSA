\name{JGGM}
\alias{JGGM}
\title{ Joint estimation of Multiple Gaussian Graphical Models}
\description{Infer networks from Multiple Gaussian data from differnt groups using our proposed algorithm.}

\usage{
JGGM(data,ALPHA1=0.05,ALPHA2=0.01)
}

\arguments{
  \item{ data }{ a list of \emph{n}x\emph{p} data matrices. n can be different for each dataset but p should be the same. }
  \item{ ALPHA1 }{ The significance level of correlation screening. In general, a high significance level of correlation screening will lead to 
 a slightly large separator set S_{ij}, which reduces the risk of missing some important variables in 
 the conditioning set. Including a few false variables in the conditioning set will not hurt much the 
 accuracy of the \eqn{\psi}-partial correlation coefficient.}
  \item{ ALPHA2  }{ The significance level of \eqn{\psi} screening.}
}
\details{
This is the function that can jointly estimate multiple GGMs which can integrate the information throughtout all datasets. The method mainly consists three steps: (i) separate estimation of \eqn{\psi}-scores for each dataset, (ii) identifies possible changes of each edge across different groups and integrate the \eqn{\psi} scores across different groups simultaneously and (iii) apply multiple hypothesis test to identify edges using integrated \eqn{\psi} scores. 
}
\value{
A list of two elements:
  \item{A}{An array of multiple adjacency matrices of networks which is a \emph{M}x\emph{p}x\emph{p} array. \emph{M} is the number of dataset groups, \emph{p} is the dimension of variables in each group.}
  \item{score}{Estimated integrative \eqn{\psi} scores matrix for all pairs of different datasets. The first two columns denote the pair indices of variables i and j and the rest columns denote the Estimated integrative \eqn{\psi} scores for this pair in different groups. }
%% ...
}

\author{ Bochao Jia\email{jbc409@ufl.edu} and Faming Liang}

\examples{
#library(equSA)
#data(SR0)
#data(TR0)
#data_all <- vector("list",2)
#data_all[[1]] <- SR0
#data_all[[2]] <- TR0
#JGGM(data_all,ALPHA1=0.05,ALPHA2=0.01)
         
}

\references{
Liang, F., Song, Q. and Qiu, P. (2015). An Equivalent Measure of Partial Correlation Coefficients for High Dimensional Gaussian Graphical Models. J. Amer. Statist. Assoc., 110, 1248-1265.

Liang, F. and Zhang, J. (2008) Estimating FDR under general dependence using stochastic approximation.  Biometrika,  95(4), 961-977.

Jia, B. and Liang, F. (2017) Joint Estimation of Multiple Gaussian Graphical Models via Multiple Hypothesis Tests (preparing)

}

\keyword{JGGM}

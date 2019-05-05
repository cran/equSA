\name{JMGM}
\alias{JMGM}
\title{ Joint Mixed Graphical Models}
\description{Infer network structures from multiple datasets with mixed types of variables and edge restrictions option available.}

\usage{
JMGM(data,ALPHA1=0.05,ALPHA2=0.01,restrict=FALSE,parallel=FALSE,nCPUs)
}

\arguments{
  \item{ data }{ a list of \eqn{n}x\eqn{p} data matrices. \eqn{n} can be different for each dataset but \eqn{p} should be the same. }
  \item{ ALPHA1 }{ The significance level of correlation screening. In general, a high significance level of correlation screening will lead to 
 a slightly large separator set \eqn{S_{ij}}, which reduces the risk of missing some important variables in 
 the conditioning set. Including a few false variables in the conditioning set will not hurt much the 
 accuracy of the \eqn{\psi}-partial correlation coefficient.}
  \item{ ALPHA2  }{ The significance level of \eqn{\psi} screening.}
  \item{restrict}{Should edge restriction applied? (logical). If \code{TRUE}, we assume that there should be no edge among binary variables. The default is \code{FALSE}.}
  \item{parallel}{Should parallelization be used? (logical), default is \code{FALSE}.}
  \item{nCPUs}{Number of cores used for parallelization. Recommend to be equal to the number of datasets.}
}
\details{
This is the function that can jointly estimate multiple graphical models with mixed types of data and also consider the edge restriction scenarios. The method has three novelties: First, the proposed method resolves the conditional independence information using a \eqn{p}-learning algorithm and therefore can be applied to the mixed types of random variables. Second, the proposed method can construct networks with restricted edges determined by some preliminary knowledges. Third, the proposed method involves a Fast Bayesian joint estimation method which works on edge-wise scores and can achieve both fast and accurate integration performance for constructing multiple networks. See Jia and Liang (2018).
}
\value{
A list of three elements:
  \item{A}{An array of multiple adjacency matrices of networks which is a \eqn{M}x\eqn{p}x\eqn{p} array. \eqn{M} is the number of dataset groups, \eqn{p} is the dimension of variables in each group.}
    \item{score.sep}{Separately estimated \eqn{\psi} scores matrix for all pairs in multiple datasets. The first two columns denote the pair indices of variables \eqn{i} and \eqn{j} and the rest columns denote the estimated \eqn{\psi} scores for this pair in different groups. }
  \item{score.joint}{Estimated integrative \eqn{\psi} scores matrix for all pairs in multiple datasets. The first two columns denote the pair indices of variables \eqn{i} and \eqn{j} and the rest columns denote the estimated integrative \eqn{\psi} scores for this pair in different groups. }
%% ...
}

\author{ Bochao Jia\email{jbc409@gmail.com} and Faming Liang}

\examples{

\dontshow{
library(equSA)
data1 <- DAGsim(n=200, p=20, p.binary=0, type="AR(2)")$data
data_all <- vector("list",1)
data_all[[1]] <- data1
JMGM(data_all,ALPHA1=0.5,ALPHA2=0.9)

}
\donttest{
library(equSA)
data1 <- DAGsim(n=200, p=100, type="AR(2)")$data
data2 <- DAGsim(n=200, p=100, type="AR(2)")$data
data_all <- vector("list",2)
data_all[[1]] <- data1
data_all[[2]] <- data2
JMGM(data_all,ALPHA1=0.1,ALPHA2=0.05,parallel=TRUE,nCPUs=2)

}


         
}

\references{
Jia, B., and Liang, F. (2018) Joint Estimation of Restricted Mixed Graphical Models. manuscript.

}

\keyword{JMGM}

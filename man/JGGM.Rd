\name{JGGM}
\alias{JGGM}
\title{ Joint estimation of Multiple Gaussian Graphical Models}
\description{Infer networks from Multiple Gaussian data from differnt groups using our proposed fast Bayesian integrative method.}

\usage{
JGGM(data,ALPHA1=0.05,ALPHA2=0.01,structure = "temporal",parallel=FALSE,nCPUs)
}

\arguments{
  \item{ data }{ a list of \eqn{n}x\eqn{p} data matrices. \eqn{n} can be different for each dataset but \eqn{p} should be the same. }
  \item{ ALPHA1 }{ The significance level of correlation screening. In general, a high significance level of correlation screening will lead to 
 a slightly large separator set \eqn{S_{ij}}, which reduces the risk of missing some important variables in 
 the conditioning set. Including a few false variables in the conditioning set will not hurt much the 
 accuracy of the \eqn{\psi}-partial correlation coefficient.}
  \item{ ALPHA2  }{ The significance level of \eqn{\psi} screening.}
    \item{ structure  }{ The depedent structure of multiple networks, either "temporal" or "spatial". The default is "temporal".}
  \item{parallel}{Should parallelization be used? (logical), default is \code{FALSE}.}
  \item{nCPUs}{Number of cores used for parallelization. Recommend to be equal to the number of datasets.}
}
\details{
This is the function that can jointly estimate multiple GGMs which can integrate the information throughtout all datasets. The method mainly consists three steps: (i) separate estimation of \eqn{\psi}-scores for each dataset, (ii) identifies possible changes of each edge across different groups and integrate the \eqn{\psi} scores across different groups simultaneously and (iii) apply multiple hypothesis test to identify edges using integrated \eqn{\psi} scores. See Jia, B., et al (2018).
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
data1 <- GauSim(100,15)$data
data2 <- GauSim(100,15)$data
data_all <- vector("list",2)
data_all[[1]] <- data1
data_all[[2]] <- data2
JGGM(data_all,ALPHA1=0.5,ALPHA2=0.5)

}
\donttest{
library(equSA)
data(SR0)
data(TR0)
data_all <- vector("list",2)
data_all[[1]] <- SR0
data_all[[2]] <- TR0
JGGM(data_all,ALPHA1=0.05,ALPHA2=0.05)

}


         
}

\references{
Jia, B., and Liang, F. (2018). A Fast Hybrid Bayesian Integrative Learning of Multiple Gene Regulatory Networks for Type 1 Diabetes. Submitted to Biostatistics.

}

\keyword{JGGM}

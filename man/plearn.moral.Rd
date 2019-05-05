\name{plearn.moral}
\alias{plearn.moral}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Learning Moral graph based on \eqn{p}-learning algorithm.
}
\description{
Construct moral graph of Bayeisan network for mixed types of random varaibles based on \eqn{p}-learning algorithm. Each variable in the dataset can be either binary or Gaussian distributed.
}
\usage{
plearn.moral(data, gaussian.index = NULL, binary.index = NULL,
alpha1 = 0.1, alpha2 = 0.02, restrict = FALSE, score.only=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{The data matrix, of dimensions \eqn{n}x\eqn{p}. Each row is an observation vector and each column is a variable.
}
  \item{gaussian.index}{The index vector of Gaussian nodes. The default value is \code{NULL}. If not specified, the system will automatically determine the index for each variable. }
  \item{binary.index}{The index vector of binary nodes. The default value is \code{NULL}. If not specified, the system will automatically determine the index for each variable.}
  \item{alpha1}{The significant level of correlation screening in \eqn{p}-learning algorithm. The default value is 0.1.}
  \item{alpha2}{The significant level of partial correlation screening in \eqn{p}-learning algorithm. The dafault value is 0.02.}
  \item{restrict}{Should edge restriction applied? (logical). If \code{TRUE}, we assume that there should be no edge among binary variables. The default is \code{FALSE}.}
  \item{score.only}{ If \code{TRUE}, it only reports \eqn{z}-scores for all pair of variables. The default is \code{FALSE}.}
}
\details{
This is the function that implements the \eqn{p}-learning algorithm for learning moral graph of Bayesian Network with mixed type of random variables. 
}
\value{
A list of two objects. 
  \item{moral.matrix}{The estimated adjacency matrix of moral graph.}
  \item{score}{The estimated \eqn{z}-scores for all pair of variables.}
}
\references{

Xu, S., Jia, B., and Liang, F. (2018). Learning Moral Graphs in Construction of High-Dimensional Bayesian Networks for Mixed Data. Submitted. 
  
}
\author{
Suwa Xu Bochao Jia and Faming Liang
}

\examples{
\dontshow{
 library(equSA)
 data.graph <- DAGsim(n = 500, p = 15, p.binary=0, type="AR(2)")$data
 plearn.moral(data.graph, alpha1 = 0.5, alpha2 = 0.5)
}
\donttest{
 library(equSA)
 data.graph <- DAGsim(n = 200, p = 100, type="AR(2)", p.binary = 50)$data
 plearn.moral(data.graph, alpha1 = 0.1, alpha2 = 0.02)
}




}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ plearn.moral }


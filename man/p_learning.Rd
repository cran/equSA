\name{p_learning}
\alias{p_learning}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construct Bayesian Network based on p-learning algorithm.
}
\description{
Construct Bayeisan network for general types of random varaibles based on \emph{p}-learning algorithm.
}
\usage{
p_learning(data, gaussian.index, binary.index, poisson.index, 
alpha1 = 0.1, alpha2 = 0.02, alpha3 = 0.02)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{The data matrix, of dimensions \emph{n}x\emph{p}. Each row is an observation vector. 
}
  \item{gaussian.index}{The index vector of continuous nodes. The default value is \code{NULL}.}
  \item{binary.index}{The index vector of binary nodes. The default value is \code{NULL}.}
  \item{poisson.index}{The index vector of poisson nodes. The default value is \code{NULL}.}
  \item{alpha1}{The significant level of step(a) of \emph{p}-screening method. The default value is 0.1.}
  \item{alpha2}{The significant level of step(c) of \emph{p}-screening method. The dafault value is 0.02.}
  \item{alpha3}{The significant level of solving Markov Blankets. The default value is 0.02.}
  
}
\details{
This is the function that implements the p-learning algorithm. 
}
\value{
A list of one object. 
  \item{PDAG}{The derived partial directed acyclic graph.}
}
\references{

Suwa, Xu and Faming, Liang (2017). Learning High-Dimensional Bayesian Networks for General Types of
Random Variables. Submitted to Biometrika. 
}
\author{
Suwa Xu and Faming Liang
}

\examples{
 #library(equSA)
 #data(mixed3000)
 #pdag3000 <- p_learning(data =mixed3000$data, gaussian.index = 
 #mixed3000$gaussian.index,binary.index <- mixed3000$binary.index)$PDAG



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ p_learning }


\name{DAGsim}
\alias{DAGsim}
\title{Simulate a directed acyclic graph with mixed data (gaussian and binary)}
\description{Simulate a directed acyclic graph with mixed data (gaussian and binary).}

\usage{
DAGsim(n, p, sparsity = 0.02,  p.binary, type="AR(2)", verbose = TRUE)
}

\arguments{
  \item{n}{Number of observations.}
  \item{p}{Number of variables. Not applicable to the graph of \code{"alarm"} type.}
\item{sparsity}{Sparsity of the graph in the \code{"random"} type, the default value is 0.02. Not applicable to other types.}
\item{p.binary}{Number of binary variables. Not applicable to the graph of \code{"alarm"} type. The default value is \code{p/2}.}

\item{type}{The graph structure with 3 options: \code{"random"}, \code{"alarm"} and \code{"AR(2)"} (default).}

\item{verbose}{If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.}
}

\details{
Given the type of graph, the patterns are generated as below:\cr\cr
(I) \code{"random"}: Each pair of off-diagonal elements are randomly set \code{edgematrix[i,j]=1} for \code{i < j} with probability \code{sparsity}, and \code{0} otherwise. It results in about \code{p*(p-1)*sparsity/2} edges in the graph.\cr\cr
(II)\code{"AR(2)"}: The off-diagonal elements are set to be \code{theta[i,j]=1} if \code{i<j} and \code{|i-j|<=2} and \code{0} otherwise.\cr\cr
(III) \code{"alarm"}: The graph structure is directly borrowed from package \pkg{'bnlearn'}, which has 37 variables with 46 edges. See \pkg{'bnlearn'} for more detail.
}

\value{
A list of five objects. 
\item{edgematrix}{A \eqn{p}x\eqn{p} matrix which indicates the true structure of directed acyclic graph. If the (i,j)th element is equal to 1, there exists a directed edge from \eqn{X_i} to \eqn{X_j}.}
\item{data}{The simulated dataset in a \eqn{n}x\eqn{p} matrix.}
\item{moral.matrix}{The simulated adjacency matrix of the moral graph, which is the undircted version of Bayesian network.}
\item{gaussian.index}{The index of Gaussian variables.}
\item{binary.index}{The index of binary variables.}
}
\references{

Kalisch, M., and Buhlmann, P. (2007). Estimating high-dimensional directed acyclic graphs with the PC-algorithm. Journal of Machine Learning Research, 8(Mar), 613-636.

Xu, S., Jia, B., and Liang, F. (2018). Learning Moral Graphs in Construction of High-Dimensional Bayesian Networks for Mixed Data. Submitted. 
  

I. A. Beinlich, H. J. Suermondt, R. M. Chavez, and G. F. Cooper. The ALARM Monitoring System: A Case Study with Two Probabilistic Inference Techniques for Belief Networks. In Proceedings of the 2nd European Conference on Artificial Intelligence in Medicine, pages 247-256. Springer-Verlag, 1989.

}
\author{
Suwa Xu, Bochao Jia and Faming Liang
}

\examples{
\dontshow{
library(equSA)
DAGsim(n=300, type="alarm")
}
\donttest{
library(equSA)
DAGsim(n=300, p=100, type="AR(2)", p.binary=50)
}

}
\keyword{DAGSim}


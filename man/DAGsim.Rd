\name{DAGsim}
\alias{DAGsim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Simulate a directed acyclic graph with mixed data (continuous and binary)
}
\description{
Simulate a directed acyclic graph with mixed data (continuous and binary).
}
\usage{
DAGsim(n, p, sparsity = 0.02, p.binary)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{x} here~~
number of observations.
}
  \item{p}{
%%     ~~Describe \code{x} here~~
number of variables.
}

  \item{sparsity}{
%%     ~~Describe \code{x} here~~
sparsity of the graph.
}

  \item{p.binary}{
%%     ~~Describe \code{x} here~~
number of binary variables.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
The default value of sparsity is 0.02.
}
\value{


A list of four objects. 
\item{Adjacency.matrix}{\emph{p}x\emph{p} The simulated adjacency matrix which indicates the true structure of directed acyclic graph. If the (i,j)th element is equal to 1, there exists a directed edge from \eqn{X_i} to eqn{X_j}.}
\item{Data}{The simulated dataset in a \emph{n}x\emph{p} matrix.}

\item{gaussian.index}{The index of continuous variables.}

\item{binary.index}{The index of binary variables.}
}
\references{
%% ~put references to the literature/web site here ~
Kalisch, Markus and B{\"u}hlmann, Peter (2007). Estimating high-dimensional directed acyclic graphs with the PC-algorithm. Journal of Machine Learning Research, 8, 613-636.

Suwa, Xu and Faming, Liang (2017). Learning High-Dimensional Bayesian Networks for General Types of
Random Variables. Submitted to Biometrika. 


}
\author{
%%  ~~who you are~~
Suwa Xu, Faming Liang
}

\examples{
# library(equSA)
# set.seed(3)
# dagsim <- DAGsim(n = 3000, p = 100, sparsity = 0.02, p.binary = 50)
# data3000 <- dagsim$data
# cont_index  <- dagsim$gaussian.index
# binary_index<- dagsim$binary.index
# truegraph <- dagsim$Adjacency.matrix
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{DAGSim}


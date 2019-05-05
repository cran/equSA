\name{ContSim}
\alias{ContSim}
\title{ A simulation method for generating count data from multivariate Zero-Inflated Negative Binomial distributions}
\description{Implements the data generation from multivariate Zero-Inflated Negative Binomial (ZINB) distributions with different graph structures, including \code{"random"}, \code{"hub"}, \code{"cluster"}, \code{"AR(2)"} and \code{"scale-free"}.}

\usage{
ContSim(n, p,  v = NULL, u = NULL, g = NULL, prob = NULL, vis = FALSE, verbose = TRUE, 
graph.type="AR(2)", k=3.30, lambda=515, omega=0.003,lower.tail = TRUE, log.p = FALSE)
}

\arguments{
  \item{ n }{ The number of observations (sample size).}
  \item{ p }{ The number of variables (dimension).}
   \item{graph.type}{
The graph structure with 4 options: \code{"random"}, \code{"hub"}, \code{"cluster"}, \code{"AR(2)"} and \code{"scale-free"}.
}
  \item{v}{
The off-diagonal elements of the precision matrix, controlling the magnitude of partial correlations with \code{u}. The default value is 0.3.
}
  \item{u}{
A positive number being added to the diagonal elements of the precision matrix, to control the magnitude of partial correlations. The default value is 0.1.
}
  \item{g}{
For \code{"cluster"} or \code{"hub"} graph, \code{g} is the number of hubs or clusters in the graph. The default value is about \code{d/20} if \code{d >= 40} and 2 if \code{d < 40}. NOT applicable to \code{"random"} and \code{"AR(2)"} graph.
}
  \item{prob}{
For \code{"random"} graph, it is the probability that a pair of nodes has an edge. The default value is \code{3/d}. For \code{"cluster"} graph, it is the probability that a pair of nodes has an edge in each cluster. The default value is \code{6*g/d} if \code{d/g <= 30} and 0.3 if \code{d/g > 30}. NOT applicable to \code{"hub"} or \code{"AR(2)"} graphs.
}
  \item{vis}{
Visualize the adjacency matrix of the true graph structure, the graph pattern, the covariance matrix and the empirical covariance matrix. The default value is \code{FALSE}
}
  \item{verbose}{
If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
}
 
  \item{k}{dispersion parameter of ZINB distribution, default of 3.30.}

  \item{lambda}{vector of (non-negative) means of ZINB distribution, default of 515.}

  \item{omega}{zero-inflation parameter of ZINB distribution, default of 0.003.}

  \item{lower.tail}{logical; if TRUE (default), probabilities are \code{P[X <= x]}, otherwise, \code{P[X> x]}.}
  
  \item{log.p}{logical; if TRUE, probabilities \code{p} are given as \code{log(p)}.}
}
\details{
This is the function that can generate dataset from multivariate Zero-Inflated Negative Binomial distributions with different graph structures, including \code{"random"}, \code{"hub"}, \code{"cluster"}, \code{"AR(2)"} and \code{"scale-free"}.

Given the adjacency matrix \code{theta}, the graph patterns are generated as below:\cr\cr
(I) \code{random}: Each pair of off-diagonal elements are randomly set \code{theta[i,j]=theta[j,i]=1} for \code{i!=j} with probability \code{prob}, and \code{0} other wise. It results in about \code{d*(d-1)*prob/2} edges in the graph.\cr\cr
(II)\code{hub}:The row/columns are evenly partitioned into \code{g} disjoint groups. Each group is associated with a "center" row \code{i} in that group. Each pair of off-diagonal elements are set \code{theta[i,j]=theta[j,i]=1} for \code{i!=j} if \code{j} also belongs to the same group as \code{i} and \code{0} otherwise. It results in \code{d - g} edges in the graph.\cr\cr
(III)\code{cluster}:The row/columns are evenly partitioned into \code{g} disjoint groups. Each pair of off-diagonal elements are set \code{theta[i,j]=theta[j,i]=1} for \code{i!=j} with the probability \code{prob}if both \code{i} and \code{j} belong to the same group, and \code{0} other wise. It results in about \code{g*(d/g)*(d/g-1)*prob/2} edges in the graph.\cr\cr
(IV)\code{AR(2)}: The off-diagonal elements are set to be \code{theta[i,j]=0.5} if \code{|i-j|=1}, \code{theta[i,j]=0.05} if \code{|i-j|=2} and \code{0} other wise.\cr\cr
(V) \code{scale-free}: The graph is generated using B-A algorithm. The initial graph has two connected nodes and each new node is connected to only one node in the existing graph with the probability proportional to the degree of the each node in the existing graph. It results in \code{d} edges in the graph.

The adjacency matrix \code{theta} has all diagonal elements equal to \code{0}. To obtain a positive definite precision matrix, the smallest eigenvalue of \code{theta*v} (denoted by \code{e}) is computed. Then we set the precision matrix equal to \code{theta*v+(|e|+0.1+u)I}. The covariance matrix is then computed for generating multivariate ZINB dataset.

The default values for parameters \code{k}, \code{lambda} and \code{omega} of ZINB distribution are estimated from a real TCGA dataset. See Jia.B et al(2017) for more detail.
}
\value{
A list of two elements:
  \item{data}{The simulated count dataset in a \eqn{n}x\eqn{p} matrix.}
  \item{Adj}{\eqn{p}x\eqn{p} The adjacency matrix of true graph structure (in sparse matrix representation) for the generated data}
%% ...
}

\author{ Bochao Jia\email{jbc409@gmail.com}}

\examples{
library(equSA)
ContSim(100,200)

}

\references{

Jia, B., Xu, S., Xiao, G., Lamba, V., Liang, F. (2017) Inference of Genetic Networks from Next Generation Sequencing Data. Biometrics.

T. Zhao and H. Liu.(2012) The huge Package for High-dimensional Undirected Graph Estimation in R. Journal of Machine Learning Research.
 
Yahav, I., and Shmueli, G. (2012). On generating multivariate Poisson data in management science applications. Applied Stochastic Models in Business and Industry, 28(1), 91-102. 

}

\keyword{ContSim}

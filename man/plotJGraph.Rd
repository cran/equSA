\name{plotJGraph}
\alias{plotJGraph}
\title{Plot Networks}
\description{Plot multiple networks with specific layout.}
\usage{
plotJGraph(A,fn="Net",th = 1e-06, mylayout = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{An array of multiple adjacency matrices of networks to be plotted which is a \eqn{M}x\eqn{p}x\eqn{p} array. \eqn{M} is the number of dataset groups, \eqn{p} is the dimension of variables in each group.}
  \item{fn}{file name to save the network plots. Default to be an string called "Net". NOTE: It should be file name for PDF file.}
  \item{th}{numeric value, default to 1e-06. To specify the threshold if the estimated coefficient between two variables is to be considered connected.}
  \item{mylayout}{graph layout to draw networks, default to \code{NULL}.}
}
\details{This function serves as the alternative plotting function to allow users to plot multiple networks with specific layout, such as plotting the simulated networks.}
\value{
Returns the multiple layout objects from \pkg{igraph} package - numeric matrix of two columns and the rows with the same number as the number of vertices.}

\author{ Bochao Jia\email{jbc409@gmail.com}, Faming liang}

\examples{
\dontshow{
library(equSA)
A1 <- GauSim(100,200,graph="scale-free")$theta
A2 <- GauSim(100,200,graph="random")$theta
A <- array(0,dim = c(2,200,200))
plotJGraph(A)
}
\donttest{

library(equSA)
data(SR0)
data(TR0)
data_all <- vector("list",2)
data_all[[1]] <- SR0
data_all[[2]] <- TR0
A <- JGGM(data_all,ALPHA1=0.05,ALPHA2=0.01)$Array
plotJGraph(A)
}

}


\references{
Jia, B., and Liang, F. (2018). Learning Multiple Gene Regulatory Networks in Type 1 Diabetes through a Fast Bayesian Integrative Method. Submitted to Journal of Statistical Computing.

}

\keyword{plotJGraph}
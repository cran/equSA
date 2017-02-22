\name{plotJGraph}
\alias{plotJGraph}
\title{Plot Networks}
\description{Plot multiple networks with specific layout.}
\usage{
plotJGraph(A,fn="Net",th = 1e-06, mylayout = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{An array of multiple adjacency matrices of networks to be plotted which is a \emph{M}x\emph{p}x\emph{p} array. \emph{M} is the number of dataset groups, \emph{p} is the dimension of variables in each group.}
  \item{fn}{file name to save the network plots. Default to be an string called "Net". NOTE: It should be file name for PDF file.}
  \item{th}{numeric value, default to 1e-06. To specify the threshold if the estimated coefficient between two variables is to be considered connected.}
  \item{mylayout}{graph layout to draw networks, default to \code{NULL}.}
}
\details{This function serves as the alternative plotting function to allow users to plot multiple networks with specific layout, such as plotting the simulated networks.}
\value{
Returns the multiple layout objects from \code{igraph} package - numeric matrix of two columns and the rows with the same number as the number of vertices.}

\author{ Bochao Jia, Faming liang\email{faliang@ufl.edu}}

\examples{
#library(equSA)
#data(SR0)
#data(TR0)
#data_all <- vector("list",2)
#data_all[[1]] <- SR0
#data_all[[2]] <- TR0
#A <- JGGM(data_all,ALPHA1=0.05,ALPHA2=0.01)$Array
#plotJGraph(A)
}


\references{
Allen, G.I., and Liu, Z. (2012). A Log-Linear graphical model for inferring genetic networks from high-throughput sequencing data. \emph{The IEEE International Conference on Bioinformatics and Biomedicine (BIBM 2012)}.

Jia, B. and Liang, F. (2017) Joint Estimation of Multiple Gaussian Graphical Models via Multiple Hypothesis Tests (preparing)

}

\keyword{plotJGraph}
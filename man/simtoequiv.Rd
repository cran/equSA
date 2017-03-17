\name{simtoequiv}
\alias{simtoequiv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Transform a directed acyclic graph into an equivalent correct graph.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A correct graph is specified  by its adjacencies and V-structures only.
}

\usage{
simtoequiv(edgematrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{edgematrix}{The simulated true graph.

}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~

%}
\value{
A list of one object. 
 \item{PDAG }{The equivalent correct graph.}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{

Jean-Philippe, Pellet and Andr{\'e},Elisseeff (2008). Using Markov blankets for causal structure learning. Journal of Machine Learning Research, 9, 1295-1342. 

Suwa, Xu and Faming, Liang (2017). Learning High-Dimensional Bayesian Networks for General Types of
Random Variables. Submitted to Biometrika. 

}
\author{
Suwa Xu, Faming Liang
}

\examples{
# library(equSA)
# load("mixed3000.rda")
# equiv_graph <- simtoequiv(mixed3000$Adjacency.matrix)$PDAG

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ simtoequiv}

\name{Mulpval}
\alias{Mulpval}
\title{ Multiple hypothesis tests for \eqn{p} values}
\description{Conduct multiple hypothesis tests from \eqn{p} values.}

\usage{
Mulpval(pvalue, ALPHA2=0.05,GRID=2,iteration=100)
}

\arguments{
  \item{ pvalue }{ A vector of \eqn{p} values.}
 \item{ ALPHA2  }{ The significance level of screening, default of 0.05.}
    \item{GRID}{The number of components for the \eqn{z}-scores. The default value is 2.}
\item{iteration}{Number of iterations for screening. The default value is 100.}
}
\details{
This is the function that conduct multiple hypothesis test for \eqn{p} values.
}
\value{
  \item{qqqscore}{The threshold of \eqn{p} value which indicates that \eqn{p} values are not larger than the threshold are considered significance and larger otherwise.
%% ...
}
}

\author{ 
Bochao Jia, Faming liang\email{fmliang@purdue.edu}
}

\examples{
library(equSA)
pvalue <- c(runif(20,0,0.001),runif(200,0,1))
Mulpval(pvalue,ALPHA2=0.05)
}

\references{

Liang, F. and Zhang, J. (2008) Estimating FDR under general dependence using stochastic approximation.  Biometrika,  95(4), 961-977.
}

\keyword{Mulpval}

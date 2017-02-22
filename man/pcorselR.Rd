\name{pcorselR}
\alias{pcorselR}
\title{ Multiple hypothesis test}
\description{Infer networks from \eqn{\psi} scores using multiple hypothesis test in \eqn{\psi} screening procedure.}

\usage{
pcorselR(score, ALPHA2=0.05)
}

\arguments{
  \item{ score }{ \eqn{\psi} score matrix which has 3 columns. The first two columns denote the pair of variables i and j and the last column denote the calculated \eqn{\psi} scores for this pair.}
 \item{ ALPHA2  }{ The significance level of \eqn{\psi} screening, default of 0.05.}
}
\details{
This is the function that conduct multiple hypothesis test for \eqn{\psi} scores, thus we called it \eqn{\psi} screening procedure.
}
\value{
  \item{qqqscore}{The threshold value of \eqn{\psi} scores which indicates that if one pair of variables has larger \eqn{\psi} scores than this threshold value in the \eqn{\psi} score matrix, this pair is considered as connected, i.e there is an edge between this pair of variables.
%% ...
}
}

\author{ 
Bochao Jia, Faming liang\email{faliang@ufl.edu}
}

\examples{
# library(equSA)
# data(SR0)
# U <- psical(SR0, ALPHA1=0.05)
# ##   probit transformation for psi scores ###
# z<-U[,3]
# q<-pnorm(-abs(z), log.p=TRUE)
# q<-q+log(2.0)
# s<-qnorm(q,log.p=TRUE)
# s<-(-1)*s
# U<-cbind(U[,1:2],s)
# ## subsampling for psi scores ###
# N <- length(U[,1])
# ratio<-ceiling(N/100000)
# U<-U[order(U[,3]), 1:3]
# m<-floor(N/ratio)
# m0<-N-m*ratio
# s<-sample.int(ratio,m,replace=TRUE)
# for(i in 1:length(s)) s[i]<-s[i]+(i-1)*ratio
# if(m0>0){
#   s0<-sample.int(m0,1)+length(s)*ratio
#   s<-c(s,s0)
# }
# Us<-U[s,]
# y <- round(Us,6)
# ##  multiple hypothesis tests ###
# pcorselR(y,ALPHA2=0.05)       
}

\references{
Liang, F., Song, Q. and Qiu, P. (2015). An Equivalent Measure of Partial Correlation Coefficients for High Dimensional Gaussian Graphical Models. J. Amer. Statist. Assoc., 110, 1248-1265.

Liang, F. and Zhang, J. (2008) Estimating FDR under general dependence using stochastic approximation.  Biometrika,  95(4), 961-977.
}

\keyword{pcorselR}

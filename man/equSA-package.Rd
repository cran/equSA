\name{equSA-package}
\alias{equSA-package}
\docType{package}
\title{

Graphical model has been widely used in may scientific fileds to describe the conditional independent relationships for a large set of random variables. Through this package, we provide tools to learn both undirected graph (Markov Random Field) and directed acyclic graph (Bayesian Network).
p
}
\description{
The package contains two parts, learning undirected graph and directed acyclic graph. 

In the first part, the package provides an equivalent measure of partial correlation coefficients for high-dimensional Gaussian Graphical Models to learn and visualize the underlying relationships between variables from single or multiple datasets. The package also provides the method for constructing networks for Next Generation Sequencing Data. Besides, it includes the method for jointly estimating Gaussian Graphical Models of multiple datasets.

In the second part, the package implements the p-learning algorithm which is used to learn Bayesian networks for general types of random variables. 
}
\details{
\tabular{rl}{
Package: \tab equSA\cr
Type: \tab Package\cr
Version: \tab 1.1.4\cr
Date: \tab 22018-01-09\cr
License: \tab GPL-2\cr
}

We propose an equvalent mearsure of partial correlation coeffient estimator called \eqn{\psi} estimators which enable us to estimate these networks via sparse, high-dimensional undirected graphical models. (Liang, F et al, 2015)

\bold{}

Here, we provide the community a convenient and useful tool to learn a Gaussian Graphical Models.

To estimate the network structures from Gaussian distributed data with this package, users simply need to specify the \code{"method"} in the main function, 
for example \code{equSAR(data,...)} to fit GGM to get the estimated adjacency matrix.
\bold{}

In this package, we also provide the code for combining Networks from two different dataset \code{combineR(data1,data2,...)} and the
code for detecting difference between two Networks, for example \code{diffR(data1,data2,...)}. data1 and data2 should share the same dimension of variables (p) but allow have different samples (n). 

\bold{}

This package also implement the Algorithm 17.1 of Friedman et al(2001), i.e estimate the covariance and precision matrix of the data given its structure. \code{solcov(data,struct,...)}

\bold{}
Besides estimating single GGM, we also propose a joint estimation method for multiple GGM. This is achieved by \eqn{\psi}- learning algorithm for graphical model at each time point combined with an Bayesian data integration method to estimate integrative \eqn{\psi} scores. Then multiple hypothesis tests were applied to identify the edges for each pair of variables. \code{JGGM(data,...)}.

\bold{}
If the data are not Normalized, for example, the count data, we propose a random effect model-based transformation to continuized data \code{ContTran(data,...)}, and then we transform the continuized data to Gaussian via a semiparametric transformation and then apply \eqn{\psi}- learning algorithm to reconstruct networks. The proposed method is consistent, and the resulting network satisfies
the faithfulness and global Markov properties.The most common application is to estimate Gene Regulatory Networks from Next Generation Sequencing Data (Jia, B et al, 2017)


\bold{}
For learning Bayesian network, the package currently supports for Gaussian, binary and poisson data and also mixed type of data. \code{p_learning(data,...)}. The proposed algorithm provides a feasible way to describe conditional dependence relationships. A straightforward application of the Bayesian network is selection of causal features for high-dimensional generalized linear model.



}
\author{
Bochao Jia, Faming Liang, Runmin Shi, Suwa Xu
Maintainer: Bochao Jia<jbc409@ufl.edu>
}
\references{
Friedman, J., Hastie, T., & Tibshirani, R. (2001). The elements of statistical learning (Vol. 1). Springer, Berlin: Springer series in statistics.

Liang, F., Song, Q. and Qiu, P. (2015). An Equivalent Measure of Partial Correlation Coefficients for High Dimensional Gaussian Graphical Models. J. Amer. Statist. Assoc., 110, 1248-1265.<doi:10.1080/01621459.2015.1012391>

Liang, F. and Zhang, J. (2008) Estimating FDR under general dependence using stochastic approximation.  Biometrika,  95(4), 961-977.<doi:10.1093/biomet/asn036>


Liu, H., Lafferty, J. and Wasserman, L. (2009). The Nonparanormal: Semiparametric Estimation 
 of High Dimensional Undirected Graphs. Journal of Machine Learning Research , 10, 2295-2328.
 


Jia, B., Xu, S., Xiao, G., Lamba, V., Liang, F. (2017) Inference of Genetic Networks from Next Generation Sequencing Data. Biometrics.

Jia, B. and Liang, F. (2017) Joint Estimation of Multiple Gaussian Graphical Models via Multiple Hypothesis Tests (preparing)


Jean-Philippe, Pellet and Andr{\'e},Elisseeff (2008). Using Markov blankets for causal structure learning. Journal of Machine Learning Research, 9, 1295-1342. 

Suwa, Xu and Faming, Liang (2017). Learning High-Dimensional Bayesian Networks for General Types of
Random Variables. Submitted to Journal of Machine Learning. 

Kalisch, Markus and B{\"u}hlmann, Peter (2007). Estimating high-dimensional directed acyclic graphs with the PC-algorithm. Journal of Machine Learning Research, 8, 613-636.



}
\keyword{ package }
\examples{
#library(equSA)
#data(SR0)
#equSAR(SR0)
}
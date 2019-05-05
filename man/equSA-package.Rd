\name{equSA-package}
\alias{equSA-package}
\docType{package}
\title{

Graphical model has been widely used in many scientific fileds to describe the conditional independent relationships for a large set of random variables. Through this package, we provide tools to learn structure for undirected graph (Markov Random Field) and moral graph for directed acyclic graph (Bayesian Network).
}
\description{
The package provides multiple algorithms for learning high-dimensional graphical models including both undirected graph and directed acyclic graph. 
For the undirect graph, the package provides an equivalent measure of partial correlation coefficients for high-dimensional Gaussian Graphical Models. Extended methods for inferring network structures from discretevariables are also available. Moreover, we also provide some methods for estimating graphical models from multiple datasets.

For the directed acyclic graph, the package provides the \eqn{p}-learning algorithm which is used to learn moral graphs in construction of high-dimensional Bayesian Networks for mixed data.

}
\details{
\tabular{rl}{
Package: \tab equSA\cr
Type: \tab Package\cr
Version: \tab 1.2.1\cr
Date: \tab 2019-05-04\cr
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
Besides estimating single GGM, we also propose a joint estimation method for multiple GGMs. This is achieved by \eqn{\psi}- learning algorithm for graphical model at each time point combined with an Bayesian data integration method to estimate integrative \eqn{\psi} scores. Then multiple hypothesis tests were applied to identify the edges for each pair of variables. \code{JGGM(data,...)}.

\bold{}
If the data contains mixed types of variables, such as either Guassian or binary distributed. We provide a method for learning graphical models for this mixed dataset with edge restrictions option available, see \code{plearn.struct(data,...)}. We also provide a method for jointly estimation of mixed graphical model, see \code{JMGM(data,...)}.


\bold{}
If the data are not Gaussian distributed, for example, the count data, we propose a random effect model-based transformation to continuized data \code{ContTran(data,...)}, and then we transform the continuized data to Gaussian via a semiparametric transformation and then apply \eqn{\psi}- learning algorithm to reconstruct networks. The proposed method is consistent, and the resulting network satisfies
the faithfulness and global Markov properties.The most common application is to estimate Gene Regulatory Networks from Next Generation Sequencing Data (Jia, B et al, 2017).


\bold{}
If we have the data following a distinct distribution and therefore produce the heterogeneous data.
In this case, we might still be interested in constructing a single gene regulatory network for the heterogeneous data in a fashion of data integration, see \code{GGMM(data,...)}.

\bold{}
For learning high-dimensional Gaussian Graphical Models from missing data, we provide a Imputation-Regularized Optimzation (IRO) algorithm (Liang et al, 2018). See \code{GraphIRO(data,...)}.

\bold{}
For learning moral graph and markov blanket for Bayesian network, the package currently supports for Gaussian and binary data and also mixed type of data. See \code{plearn.moral(data,...)}. The proposed algorithm provides a feasible way to describe conditional dependence relationships for the directed acyclic graph.

\bold{}
To Construct confidence intervals and assess \eqn{p}-values in high-dimensional linear and generalized linear models. See \code{MNR(x,y,...)} for detail.


}
\author{
Bochao Jia, Faming Liang, Runmin Shi, Suwa Xu
Maintainer: Bochao Jia<jbc409@gmail.com>
}
\references{
Friedman, J., Hastie, T., & Tibshirani, R. (2001). The elements of statistical learning (Vol. 1). Springer, Berlin: Springer series in statistics.

Liang, F., Song, Q. and Qiu, P. (2015). An Equivalent Measure of Partial Correlation Coefficients for High Dimensional Gaussian Graphical Models. J. Amer. Statist. Assoc., 110, 1248-1265.<doi:10.1080/01621459.2015.1012391>

Liang, F. and Zhang, J. (2008) Estimating FDR under general dependence using stochastic approximation.  Biometrika,  95(4), 961-977.<doi:10.1093/biomet/asn036>


Liang, F., Jia, B., Xue, J., Li, Q., and Luo, Y. (2018). An Imputation Regularized Optimization Algorithm for High-Dimensional Missing Data Problems and Beyond. Submitted to Journal of the Royal Statistical Society Series B.

Liu, H., Lafferty, J. and Wasserman, L. (2009). The Nonparanormal: Semiparametric Estimation of High Dimensional Undirected Graphs. Journal of Machine Learning Research , 10, 2295-2328.
 

Jia, B., Xu, S., Xiao, G., Lamba, V., Liang, F. (2017) Inference of Genetic Networks from Next Generation Sequencing Data. Biometrics.

Jia, B., and Liang, F. (2018). A Fast Hybrid Bayesian Integrative Learning of Multiple Gene Regulatory Networks for Type 1 Diabetes. Submitted to Biostatistics.


Jia, B. and Liang, F. (2018). Learning Gene Regulatory Networks with High-Dimensional Heterogeneous Data. Accept by ICSA Springer Book. 

Jean-Philippe, Pellet and Andre,Elisseeff (2008). Using Markov blankets for causal structure learning. Journal of Machine Learning Research, 9, 1295-1342. 

Xu, S., Jia, B., and Liang, F. (2019). Learning Moral Graphs in Construction of High-Dimensional Bayesian Networks for Mixed Data. Neural computation, 1-32.
  
Jia, B., and Liang, F. (2018) Joint Estimation of Restricted Mixed Graphical Models. manuscript.

Liang, F., Xue, J. and Jia, B. (2018). Markov Neighborhood Regression for High-Dimensional Inference. Submitted to J. Amer. Statist. Assoc.


}
\keyword{ package }
\examples{
\dontshow{
library(equSA)
data(TR0)
subset <- TR0[1:100,1:10]
equSAR(subset,GRID=2,iteration=100)
}
\donttest{
library(equSA)
data(TR0)
subset <- TR0
equSAR(subset)
}

}
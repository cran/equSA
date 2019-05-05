\name{Cont2Gaus}
\alias{Cont2Gaus}
\title{ A transfomation from count data into Gaussian data}
\description{To transform count data into Gaussian distributed and also keep the consistency for contructing networks.}

\usage{
Cont2Gaus(iData,total_iteration=5000,stepsize=0.05)
}

\arguments{
  \item{ iData }{ a \eqn{n}x\eqn{p} count data matrix.}
  \item{ total_iteration  }{ Total iteration number for Baysian random effect model-based transformation, default of 5000.}
  \item{ stepsize  }{The stepsize of updating parameters in transformation, default of 0.05.}
}
\details{
This is the function that transform the count data into Gaussian data which include two steps. First, we do data continuized transformation \code{ContTran(data,...)} and then we apply the semiparametric transformation (Liu, H et al, 2009) provided in \pkg{huge} packages to tranform continuized data into Gaussian distributed. 

}
\value{
  \item{Gaus}{A \eqn{n}x\eqn{p} matrix of normalized data with Gaussian distribution.}

%% ...
}

\author{ Bochao Jia\email{jbc409@gmail.com} and Faming Liang}


\examples{
      \dontshow{
            library(equSA)
            data(count)
            subset <- count[1:5,1:10]
            Cont2Gaus(subset,total_iteration=10)
      }
      \donttest{
           library(equSA)
           data(count)
           Cont2Gaus(count,total_iteration=1000)
      }
}

\references{

Jia, B., Xu, S., Xiao, G., Lamba, V., Liang, F. (2017) Inference of Genetic Networks from Next Generation Sequencing Data. Biometrics.

Liu, H., Lafferty, J. and Wasserman, L. (2009). The Nonparanormal: Semiparametric Estimation of High Dimensional Undirected Graphs. Journal of Machine Learning Research , 10, 2295-2328.
}

\keyword{Cont2Gaus}

\name{solcov}
\alias{solcov}
\title{Calculate covariance matrix and precision matrix}
\description{
Calculate the adjusted covriance matrix and precision matrix given the network structure from high dimesional dataset.
}
\usage{
solcov(data, struct, tol=10^-5)
}

\arguments{
  \item{ data }{ A \emph{n}x\emph{p} data matrix.}
  \item{ struct }{ A preacquired adjacency matrix }
  \item{ tol  }{ Tolerant value, default is 10^-5 }
  
  }

\value{
A list of two elements:
 \item{ COV }{ Adjusted covriance matrix}
 \item{ PRE }{ Precision matrix }
%% ...
}
\author{ Bochao Jia& Runmin Shi \email{jbc409@ufl.edu}}


\examples{
library(equSA)
data(SR0)
data(SR0_mat)
subSR0 <- SR0[1:10,1:10]
subSR0_mat <- SR0_mat[1:10,1:10]
solcov(subSR0, subSR0_mat) 
      \dontshow{
            library(equSA)
            data(SR0)
            data(SR0_mat)
            subSR0 <- SR0[1:10,1:10]
            subSR0_mat <- SR0_mat[1:10,1:10]
            solcov(subSR0, subSR0_mat) 
      }
      \donttest{
           # library(equSA)
           # data(SR0)
           #  data(SR0_mat)
           # solcov(SR0, SR0_mat)
      }
}

\references{
Friedman, J., Hastie, T., & Tibshirani, R. (2001). The elements of statistical learning (Vol. 1). Springer, Berlin: Springer series in statistics.
}

\keyword{solcov}
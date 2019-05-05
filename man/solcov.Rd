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
  \item{ data }{ A \eqn{n}x\eqn{p} data matrix.}
  \item{ struct }{ A preacquired adjacency matrix }
  \item{ tol  }{ Tolerant value, default is 10^-5 }
  
  }

\value{
A list of two elements:
 \item{ COV }{ Adjusted covriance matrix}
 \item{ PRE }{ Precision matrix }
%% ...
}
\author{ Bochao Jia\email{jbc409@gmail.com} & Runmin Shi }


\examples{
  \dontshow{
library(equSA)
data <- GauSim(20,10)
solcov(data$data,data$theta)
      }
      \donttest{
library(equSA)
data <- GauSim(100,200)
solcov(data$data,data$theta)
      }
}

\references{
Friedman, J., Hastie, T., and Tibshirani, R. (2001). The elements of statistical learning (Vol. 1). Springer, Berlin: Springer series in statistics.
}

\keyword{solcov}
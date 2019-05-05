\name{ContTran}
\alias{ContTran}
\title{ A data continuized transformation}
\description{Transform count data into continuous data.}

\usage{
ContTran(iData,total_iteration=5000,stepsize=0.05)
}

\arguments{
  \item{ iData }{ a \eqn{n}x\eqn{p} count data matrix.}
  \item{ total_iteration  }{ total iteration number for Baysian random effect model-based transformation, default of 5000.}
  \item{ stepsize  }{The stepsize of updating parameters in transformation, default of 0.05.}
}
\details{
This is the function that transform the count data into continuized data.
}
\value{
  \item{continuz}{\eqn{n}x\eqn{p} matrix of continuized data.}

%% ...
}

\author{ Bochao Jia\email{jbc409@gmail.com}, Suwa Xu and Faming Liang}

\examples{
      \dontshow{
            library(equSA)
            data(count)
            subset <- count[1:5,1:10]
            ContTran(subset,total_iteration=10)
      }
      \donttest{
           library(equSA)
           data(count)
           ContTran(count,total_iteration=1000)
      }
}

\references{

Jia, B., Xu, S., Xiao, G., Lamba, V., Liang, F. (2017) Inference of Genetic Networks from Next Generation Sequencing Data. Biometrics.
}

\keyword{ContTran}

\name{ContTran}
\alias{ContTran}
\title{ A data continuized transformation}
\description{To transform count data into continuous data.}

\usage{
ContTran(iData,total_iteration=5000,stepsize=0.05)
}

\arguments{
  \item{ iData }{ a \emph{n}x\emph{p} count data matrix.}
  \item{ total_iteration  }{ total iteration number for Baysian random effect model-based transformation, default of 5000.}
  \item{ stepsize  }{The stepsize of updating parameters in transformation, default of 0.05.}
}
\details{
This is the function that transform the count data into continuized data.
}
\value{
  \item{continuz}{\emph{n}x\emph{p} matrix of continuized data.}

%% ...
}

\author{ Bochao Jia\email{jbc409@ufl.edu}, Suwa Xu and Faming Liang}

\examples{
#library(equSA)
#data(count)
#ContTran(count,total_iteration=1000)
         
}

\references{

Jia, B., Xu, S., Xiao, G., Lamba, V., Liang, F. (2017) Inference of Genetic Networks from Next Generation Sequencing Data. Biometrics, in press.
}

\keyword{ContTran}

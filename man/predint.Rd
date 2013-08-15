\name{predint}
\alias{predint}
\alias{nparpredint}
\title{Prediction intervals to contain at least k out of m future observations} 
\description{Calculates parametric and non-parametric prediction
  intervals to contain at least k out of m future observations given a observed sample}

\usage{
predint(x, k, m, level=0.95, alternative="two.sided", quantile=NULL, absError=0.001, interval=c(0, 100))
nparpredint(x, k, m, level=0.95, alternative="two.sided")
}

\arguments{
  \item{x}{A numeric vector with a sample of observations}
  \item{k}{Number of future observations contained in the interval}
  \item{m}{Size of future sample}
  \item{level}{Confidence level}
  \item{alternative}{One of "two.sided", "less", or "greater" to
    calculate two- or one-sided intervals}
  \item{quantile}{Supply a user-defined quantile for parametric intervals}
  \item{absError}{The maximum absolute error tolerated when calculating the quantile}
  \item{interval}{Root finding interval when calculating the quantile}
}


\value{ An object inheriting from class
  \code{\linkS4class{PredInterval}}, or \code{\linkS4class{nparPredInterval}}}

\seealso{\code{\link{precint}}}

\references{
  Odeh, RE (1990): 2-Sided prediction intervals to contain at least k
  out of m future observations from a normal distribution. Technometric
  32(2): 203-216.

  Danziger, L. and Davis, S.A. (1964): Tables of distribution-free tolerance limits. Annals of Mathematical Statistics 35(3):1361-1365.
}
\examples{
x <- rnorm(50)
predint(x, k=7, m=10)
nparpredint(x, k=7, m=10)
}

\keyword{htest}



\name{precint}
\alias{precint}
\alias{nparprecint}
\title{Prediction intervals to contain the mean or median of m future observations} 
\description{Calculates parametric and non-parametric prediction
  intervals to contain the mean of m future observations given a observed sample}

\usage{
precint(x, m, level=0.95, alternative="two.sided")
nparprecint(x, m, level=0.95, alternative="two.sided")
}

\arguments{
  \item{x}{A numeric vector with a sample of observations}
  \item{m}{Size of future sample}
  \item{level}{Confidence level}
  \item{alternative}{One of "two.sided", "less", or "greater" to
    calculate two- or one-sided intervals}
}



\value{ An object inheriting from class
  \code{\linkS4class{PrecInterval}}, or \code{\linkS4class{nparPrecInterval}}}

\seealso{\code{\link{predint}}}

\references{
Hahn, GJ and Meeker, WQ (1991): Statistical Intervals. Wiley, New York.

Chakraborti, S., Van der Laan, P., Van de Wiel, M.A. (2004): A class of distribution-free control charts. Applied Statistics 53(3):443-462.
}

\examples{
x <- rnorm(50)
precint(x, m=10)
nparprecint(x, m=10)
}

\keyword{htest}



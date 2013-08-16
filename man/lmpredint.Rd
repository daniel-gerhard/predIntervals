\name{lmpredint}
\alias{lmpredint}
\title{Prediction intervals to contain at least k out of m future observations for linear model fits} 
\description{Calculates parametric and non-parametric prediction
  intervals to contain at least k out of m future observations given a linear model fit of class lm}

\usage{
lmpredint(object, newdata, k, level = 0.95, alternative="two.sided", 
          quantile=NULL, absError=0.001, interval=c(0, 100))
}

\arguments{
  \item{object}{Object of class inheriting from "lm".}
  \item{newdata}{A data frame in which to look for variables with which to predict.}
  \item{k}{Number of future observations contained in the interval. k must be smaller than the number of rows of newdata.}
  \item{level}{Confidence level}
  \item{alternative}{One of "two.sided", "less", or "greater" to
    calculate two- or one-sided intervals}
  \item{quantile}{Supply a user-defined quantile}
  \item{absError}{The maximum absolute error tolerated when calculating the quantile}
  \item{interval}{Root finding interval for calculating the quantile}
}


\value{ An object inheriting from class
  \code{\linkS4class{PIlm}}}

\seealso{\code{\link{predint}}, \code{\link{predict.lm}}}

\references{
  Odeh, RE (1990): 2-Sided prediction intervals to contain at least k
  out of m future observations from a normal distribution. Technometric
  32(2): 203-216.
}

\examples{
# example from predict.lm
x <- rnorm(15)
y <- x + rnorm(15)
new <- data.frame(x = seq(-3, 3, 0.5))
m <- lm(y ~ x)

# prediction intervals to 
# contain at least 10 of 13 future observations
lmpredint(m, newdata=new, k=10)

# for a single prediction
new1 <- data.frame(x = -3)
# the prediction intervals of lmpredint
lmpredint(m, newdata=new1, k=1, absError=0.00001)
# agrees with
predict(m, newdata=new1, interval="prediction")
}

\keyword{htest}



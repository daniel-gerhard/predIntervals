\name{PInt-class}
\Rdversion{1.1}
\docType{class}
\alias{PInt-class}

\title{Class "PInt"}
\description{Objects of class "PInt" contain estimated prediction interval limits}
\section{Objects from the Class PInt}{
  Objects are created by calls of \code{\link{predint}} or \code{\link{precint}}
}
\section{Slots}{
  \describe{
    \item{\code{interval}:}{Prediction interval limits}
    \item{\code{sample}:}{Vector with observations}
    \item{\code{level}:}{Confidence level}
    \item{\code{alternative}:}{Direction of interest}
  }
}
\seealso{\code{\link{predint}}, \code{\link{precint}}}
\keyword{classes}



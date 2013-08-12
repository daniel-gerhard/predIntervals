\name{PrecInterval-class}
\Rdversion{1.1}
\docType{class}
\alias{PrecInterval-class}

\title{Class "PrecInterval"}
\description{
Objects of class "PrecInterval" contain estimated prediction interval limits
}
\section{Objects from the Class PrecInterval}{
Objects can be created by calls of \code{\link{precint}}.
}
\section{Slots}{
  \describe{
    \item{\code{m}:}{Number of future observations}
    \item{\code{interval}:}{Prediction interval limits}
    \item{\code{sample}:}{Sample observations}
    \item{\code{level}:}{Confidence level}
    \item{\code{alternative}:}{Direction of interest}
  }
}
\section{Extends}{
Class \code{"\linkS4class{PInt}"}, directly.
}
\seealso{\code{\link{precint}}}
\keyword{classes}

\name{PredInterval-class}
\Rdversion{1.1}
\docType{class}
\alias{PredInterval-class}

\title{Class "PredInterval"}
\description{
Objects of class "PredInterval" contain estimated prediction interval limits
}
\section{Objects from the Class PredInterval}{
Objects can be created by calls of \code{\link{predint}}.
}
\section{Slots}{
  \describe{
    \item{\code{m}:}{Number of future observations}
    \item{\code{k}:}{Number of future observations contained in the
      interval}
    \item{\code{interval}:}{Prediction interval limits}
    \item{\code{sample}:}{Sample observations}
    \item{\code{level}:}{Confidence level}
    \item{\code{alternative}:}{Direction of interest}
  }
}
\section{Extends}{
Class \code{"\linkS4class{PInt}"}, directly.
}
\seealso{\code{\link{predint}}}
\keyword{classes}

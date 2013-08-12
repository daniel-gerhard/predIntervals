\name{NormalPredInterval-class}
\Rdversion{1.1}
\docType{class}
\alias{NormalPredInterval-class}

\title{Class "NormalPredInterval"}
\description{
Objects of class "NormalPredInterval" contain estimated prediction interval limits
}
\section{Objects from the Class NormalPredInterval}{
Objects can be created by calls of \code{\link{predint}}.
}
\section{Slots}{
  \describe{
    \item{\code{quantile}:}{Quantile}
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
Class \code{"\linkS4class{PredInterval}"}, directly.  
Class \code{"\linkS4class{PInt}"}, by class "PredInterval", distance 2.
}
\seealso{\code{\link{predint}}}
\keyword{classes}

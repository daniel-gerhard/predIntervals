\name{NormalPrecInterval-class}
\Rdversion{1.1}
\docType{class}
\alias{NormalPrecInterval-class}

\title{Class "NormalPrecInterval"}
\description{
Objects of class "NormalPrecInterval" contain estimated prediction interval limits
}
\section{Objects from the Class NormalPrecInterval}{
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
Class \code{"\linkS4class{PrecInterval}"}, directly.  
Class \code{"\linkS4class{PInt}"}, by class "PrecInterval", distance 2.
}
\seealso{\code{\link{precint}}}
\keyword{classes}

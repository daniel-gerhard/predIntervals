\name{PIlm-class}
\Rdversion{1.1}
\docType{class}
\alias{PIlm-class}

\title{Class "PIlm"}
\description{Objects of class "PIlm" contain estimated prediction interval limits given an object of class lm}
\section{Objects from the Class PIlm}{
  Objects are created by a call of \code{\link{lmpredint}}
}
\section{Slots}{
  \describe{
    \item{\code{prediction}:}{data.frame with prediction interval limits}
    \item{\code{newdata}:}{data.frame with variables with which to predict}
    \item{\code{model}:}{object of class lm}
    \item{\code{k}:}{Number of future observations contained in the
      interval}
    \item{\code{level}:}{Confidence level}
    \item{\code{alternative}:}{Direction of interest}
  }
}
\seealso{\code{\link{lmpredint}}}
\keyword{classes}



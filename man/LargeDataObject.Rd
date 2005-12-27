\name{LargeDataObject-class}
\docType{class}
\alias{LargeDataObject-class}
\alias{show,LargeDataObject-method}
\title{Large Data Object - class}

\description{
A virtual class including the data classes \code{BeadLevelList}, \code{MAList} and \code{SegList}, all of which typically contain large quantities of numerical data in vector, matrices and data.frames.
}

\section{Methods}{
A \code{show} method is defined for objects of class \code{LargeDataObject} which uses \code{printHead} to print only the leading elements or rows of components or slots which contain large quantities of data.
}

\author{Gordon Smyth}

\examples{
#  see normalizeBetweenArrays
}

\keyword{classes}
\keyword{data}

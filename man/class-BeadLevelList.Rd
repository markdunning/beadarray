\name{BeadLevelList-class}
\docType{class}
\alias{BeadLevelList-class}
\alias{initialize,BeadLevelList-method}
\alias{[.BeadLevelList}
\alias{print.BeadLevelList}

\title{Class "BeadLevelList" }

\description{
A simple list-based class for storing red and green channel foreground and background intensities for an Illumina experiment.

}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("BeadLevelList", G, Gb, x, y, ProbeID, targets)}, but are usually created by \code{\link{readBeadLevelData}}.
}

\section{Slots/List Components}{
\code{BeadLevelList} objects can be created by \code{new("BeadLevelList",BLData)} where \code{BLData} is a list.
Objects of this class contains no slots (other than \code{.Data}), but objects should contain the following list components:
\tabular{ll}{
  \code{R}:\tab numeric matrix containing the red (foreground intensities.  Rows correspond to spots and columns to arrays.\cr
  \code{G}:\tab numeric matrix containing the green foreground intensities\cr
  \code{Rb}:\tab numeric matrix containing the red  background intensities\cr
  \code{Gb}:\tab numeric matrix containing the green background intensities
}

}
}
\section{Methods}{
  \describe{

    \item{initialize}{\code{signature(.Object = "BeadLevelList")}  }
    
  }
}

\author{Mark Dunning}
\seealso{ \code{\link{readBeadLevelData}} }
\keyword{classes}


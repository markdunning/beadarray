\name{BeadLevelList-class}
\docType{class}
\alias{BeadLevelList-class}
%\alias{BeadLevelList}
\alias{arrayNames,BeadLevelList-method}
\alias{getArrayData,BeadLevelList-method}

\title{Class "BeadLevelList" }

\description{
A class for storing red and green channel foreground and background intensities from an Illumina experiment.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("BeadLevelList")}, but are usually created by \code{\link{readIllumina}}.
}

\section{Slots/List Components}{
  Objects of this class contain the following slots
  \tabular{ll}{
  \code{beadData}:\tab an environment for storing the raw 
  bead-level data.  
  Each row correspond to a bead and columns the data.\cr
  \code{phenoData}:\tab an 'AnnotatedDataFrame' containing  
  experimental information.\cr
  \code{arrayInfo}:\tab a list containing array information.\cr
  \code{annotation}:\tab character storing annotation package 
    information.\cr
  }
}

\section{Methods}{
  \describe{
    \item{\code{arrayNames(object,arrays=NULL)}}{Returns the strip/array names from a
      \code{BeadLevelList} object for selected \code{arrays}}
    \item{\code{getArrayData(object,what="G",log=TRUE)}}{Retrieves 
      the \code{what} intensities on the \code{log} scale from the
      \code{BeadLevelList}}
  }
}
\author{Mark Dunning and Matt Ritchie}
\keyword{classes}

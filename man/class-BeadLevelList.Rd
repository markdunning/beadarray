\name{BeadLevelList-class}
\docType{class}
\alias{BeadLevelList-class}
\alias{BeadLevelList}
\alias{show,BeadLevelList-method}
\alias{dim,BeadLevelList-method}
\alias{initialize,BeadLevelList-method}
\alias{arrayNames,BeadLevelList-method}
\alias{numBeads,BeadLevelList-method}
\alias{copyBeadLevelList,BeadLevelList-method}
\alias{combineBeadLevelLists,BeadLevelList-method}
\alias{pData,BeadLevelList-method}
\alias{phenoData,BeadLevelList-method}
\alias{[[,BeadLevelList,ANY,missing-method}
\alias{getArrayData,BeadLevelList-method}
\alias{backgroundCorrect,BeadLevelList,character,double,logical-method}

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
  \code{beadAnno}:\tab a data.frame containing Bead annotation
  information.\cr
  \code{scanMetrics}:\tab a data.frame containing scanner information.
   \cr
  }
}

\section{Methods}{
  \describe{
    \item{show(BeadLevelList)}{printing method for BeadLevelList}
    \item{initialize}{\code{signature(.Object = "BeadLevelList")}}
    \item{dim}{\code{dim(object)}The dimension of the 
      \code{BeadLevelList} object}
    \item{\code{copyBeadLevelList(object)}}{Creates a new copy of a
      \code{BeadLevelList} object}
    \item{\code{arrayNames(object,arrays=NULL)}}{Returns the strip/array names from a
      \code{BeadLevelList} object for selected \code{arrays}}
    \item{\code{combineBeadLevelLists(object1,object2)}}{Combines two
      \code{BeadLevelList} objects into one}
    \item{\code{getArrayData(object,what="G",log=TRUE)}}{Retrieves 
      the \code{what} intensities on the \code{log} scale from the
      \code{BeadLevelList}}
    \item{\code{numBeads(object,arrays=NULL)}}{Returns the number of
      beads on selected \code{arrays}}
    \item{\code{pData(object)}}{Returns a data.frame with samples as
      rows, variables as columns}
    \item{\code{phenoData(object)}}{Returns an object containing
      phenotypic information on both variable values and variable
      meta-data}
  }
}
\author{Mark Dunning and Matt Ritchie}
\seealso{\code{\link{readIllumina}}}
\keyword{classes}

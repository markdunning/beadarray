\name{BeadLevelList-class}
\docType{class}
\alias{BeadLevelList-class}
\alias{coerce,BeadLevelList,exprSet2-method}
\alias{print.BeadLevelList}
\title{Red, Green Intensity List - class}

\description{
A simple list-based class for storing red and green channel foreground and background intensities for an Illumina experiment.

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


All of the matrices should have the same dimensions.
The row dimension of \code{targets} should match the column dimension of the matrices. Only the \code{G} matrix is created when single channel data is read.
}

\section{Methods}{
This class inherits directly from class \code{list} so any operation appropriate for lists will work on objects of this class.
In addition, \code{BeadLevelList} objects can be \link[limma:subsetting]{subsetted}, \link[limma:cbind]{combined} and \link[limma:merge]{merged}.
\code{BeadLevelList} objects will return dimensions and hence functions such as \code{\link[limma:dim]{dim}}, \code{\link[base:nrow]{nrow}} and \code{\link[base:nrow]{ncol}} are defined. 
\code{BeadLevelLists} also inherit a \code{\link[methods]{show}} method from the virtual class \code{\link[limma:LargeDataObject]{LargeDataObject}}, which means that \code{BeadLevelLists} will print in a compact way.

}

\author{Mark Dunning and Mike Smith}

\keyword{classes}
\keyword{data}

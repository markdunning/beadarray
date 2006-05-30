\name{cbind}
\alias{cbind.BeadLevelList}
\alias{cbind.BeadSummaryList}
\alias{rbind.BeadSummaryList}
\alias{rbind.BeadLevelList}
\title{Combine BeadLevelList or BeadSummaryList Objects}
\description{
Combine a series of \code{BeadLevelList} or combine a series of \code{BeadSummaryList}
objects. 
}
\usage{
\method{cbind}{BeadSummaryList}(\dots, deparse.level=1)
\method{rbind}{BeadSummaryList}(\dots, deparse.level=1)
}
\arguments{
  \item{\dots}{\code{BeadLevelList} objects or \code{BeadSummaryList} objects}
  \item{deparse.level}{not currently used, see \code{\link[base]{cbind}} in the base package}
}
\details{
\code{cbind} combines data objects assuming the same gene lists but different arrays.
\code{rbind} combines data objects assuming equivalent arrays, i.e., the same RNA targets, but different genes.

For \code{cbind}, the matrices of expression data from the individual objects are cbinded.
The data.frames of target information, if they exist, are rbinded.
The combined data object will preserve any additional components or attributes found in the first object to be combined.
For \code{rbind}, the matrices of expression data are rbinded while the
target information, in any, is unchanged.
}
\value{
An \code{\link[beadarray:BeadLevelList]{BeadLevelList}} or \code{\link[beadarray:BeadSummaryList]{BeadSummaryList}} object holding data from all the arrays and all genes from the individual objects.
}
\author{Gordon Smyth, modified by Mike Smith for BeadSummaryList object}
\seealso{
  \code{\link[base]{cbind}} in the base package.
  
%  \link{3.ReadingData} gives an overview of data input and manipulation functions in LIMMA.
}

\keyword{manip}

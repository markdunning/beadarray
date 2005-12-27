\name{dim}
\alias{dim.BeadLevelList}
\alias{dim.BeadSummaryList}
\alias{length.BeadLevelList}
\alias{length.BeadSummaryList}
\title{Retrieve the Dimensions of an BeadLevelList or BeadSummaryList}
\description{
Retrieve the number of rows (genes) and columns (arrays) for an BeadLevelList or BeadSummaryList.
}
\usage{
\method{dim}{BeadSummaryList}(x)
\method{length}{BeadSummaryList}(x)
}
\arguments{
  \item{x}{an object of class \code{BeadLevelList} or \code{BeadSummaryList}}
}
\details{
BeadArray data objects share many analogies with ordinary matrices in which the rows correspond to spots or genes and the columns to arrays.
These methods allow one to extract the size of microarray data objects in the same way that one would do for ordinary matrices.

A consequence is that row and column commands \code{nrow(x)}, \code{ncol(x)} and so on also work.
}
\value{
Numeric vector of length 2.
The first element is the number of rows (genes) and the second is the number of columns (arrays).
}
\author{Gordon Smyth, modified by Mike Smith for BeadSummaryList object}
\seealso{
  \code{\link[base]{dim}} in the base package.
  
}
\keyword{array}

\name{setWeights}
\alias{setWeights}
\title{Set BeadLevelList Weights}
\description{
Replaces the weights of a BeadLevelList with user-specified ones.
}
\usage{
setWeights(BLData, wts, array, combine=FALSE)
}

\arguments{
  \item{BLData}{\code{BeadLevelList}}
  \item{wts}{either a numerical vector of weights to use, or 0 or 1 to set all weights to 0 or 1.}
  \item{array}{integer specifying the strip/array to use}
  \item{combine}{logical.  If \code{TRUE}, the new weights specified by \code{wts} are combined
    with the existing weights by storing the minimum of the two for each
    bead.  If \code{FALSE} the new weights replace any existing weights.}
}

\details{
This function replaces the weights column, \code{wts}, on the specified array, with user-specified values.

Only rows with \code{wts != 1} are used in \code{createBeadSummaryData}.
}

\value{
  \code{BeadLevelList} object, with updated \code{wts} values.
}

\author{Mark Dunning}
\examples{
data(BLData)
BLData <- setWeights(BLData,1,1) ##set all weights to 1
BLData <- setWeights(BLData,0,1) ##set all weights to 0
}
\keyword{misc}

\name{setWeights}
\alias{setWeights}
\title{Set BeadLevelList Weights}
\description{
Replaces the weights of a BeadLevelList with user-specified ones.
}
\usage{
setWeights(BLData, wts, array)
}

\arguments{
  \item{BLData}{\code{BeadLevelList}}
  \item{wts}{a numerical vector }
  \item{array}{integer specifying the strip/array to use}
}

\details{
This function replaces the weights column, \code{wts}, on the specified array, with user-specified values. Only rows with \code{wts == 1} are used in CreateBeadSummary.
}

\value{
  BeadLevelList object, with updated \code{wts} values.
}

\author{Mark Dunning}
\examples{
data(BLData)
BLData <- setWeights(BLData,1,1) ##set all weights to 1
}
\keyword{misc}

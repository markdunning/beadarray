\name{numBeads}
\alias{numBeads}
\title{Gets the number of beads from a BeadLevelList object}
\description{
Retrieves the number of beads on selected strips/arrays from a \code{BeadLevelList} object.
}
\usage{
numBeads(object, arrays=NULL)
}
\arguments{
  \item{object}{\code{BeadLevelList}}
  \item{arrays}{either \code{NULL} to return the bead numbers for 
      all arrays, or a scalar or vector of integers specifying a 
      subset of strips/arrays}
}
\details{
  \code{numBeads} retrieves the number of beads on arrays from the
  arrayInfo slot.
   }
\value{
A vector containing the number of beads on individual strips/arrays.
}
\author{Matt Ritchie}
\examples{
data(BLData)
numBeads(BLData)
numBeads(BLData, arrays=2)
}
\keyword{manip}

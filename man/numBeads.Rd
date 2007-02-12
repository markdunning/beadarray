\name{numBeads}
\alias{numBeads}
\title{Gets the number of beads from a BeadLevelList object}
\description{
Retrieves the number of beads from each array in a \code{BeadLevelList} object.
}
\usage{
numBeads(object)
}
\arguments{
  \item{object}{\code{BeadLevelList}}
}
\details{
  \code{numBeads} retrieves the number of beads on each array from the
  arrayInfo slot.
   }
\value{
A vector containing the number of beads on individual arrays.
}
\author{Matt Ritchie}
\examples{
data(BLData)
numBeads(BLData)
}
\keyword{manip}

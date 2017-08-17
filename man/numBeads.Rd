\name{numBeads}
\alias{numBeads}
\title{Gets the number of beads from a beadLevelData object}
\description{
Retrieves the number of beads on selected sections from a \code{beadLevelData} object.
}
\usage{
numBeads(object, arrays=NULL)
}
\arguments{
  \item{object}{\code{beadLevelData}}
  \item{arrays}{either \code{NULL} to return the bead numbers for all arrays, or a scalar or vector of integers specifying a subset of strips/arrays}
}
\details{
  \code{numBeads} retrieves the number of beads on arrays from the arrayInfo slot.
}
\value{
A vector containing the number of beads on individual array sections.
}
\author{Matt Ritchie}
\examples{

if(require(beadarrayExampleData)){

data(exampleBLData)
numBeads(exampleBLData)
numBeads(exampleBLData, arrays=2)

}

}
\keyword{manip}
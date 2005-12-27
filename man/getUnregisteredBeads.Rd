\name{getUnregisteredBeads}
\alias{getUnregisteredBeads}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Finds unregistered beads}
\description{
Function returns a list of all beads on a particular array that were not decoded
}
\usage{
getUnregisteredBeads(BLData, array)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{an BeadLevelList object}
  \item{array}{numeric value for the array that we want to find unregistered beads on}
}
\details{
Beads which cannot be decoded are assigned a probeID of 0 and excluded from further analysis.
The function simply looks up all beads which have a probeID of 0 and returns the index of the bead (1-49777).

}
\value{
  vector of numeric specifying the indices of beads which were not decoded properly
  
}

\author{Mark Dunning}

\examples{

data(BLData)

getUnregisteredBeads(BLData, 1)

u = getUnregisteredBeads(BLData,1)

u[1:20]

BLData$probeID[33,1]

##For a quick test to see where these beads are located we can use 'plotBeadLocations'

plotBeadLocations(BLData, probeID=0, array=1, label=FALSE)

} 

\keyword{methods}
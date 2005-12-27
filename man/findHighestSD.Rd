\name{findHighestSD}
\alias{findHighestSD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Highest standard deviation.}
\description{
Function for finding which bead types on a particular array have the standard deviation of their intensities greater than some limit.
}

\usage{
findHighestSD(BLData, array, limit = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{BeadLevelList object containing bead level data}
  \item{array}{numeric value for the array that we want to use}
  \item{limit}{threshold value for the standard deviation of beads}
}
\details{
The standard deviation for bead intensities of each bead type is calculated in turn and if this value exceeds the defined limit then we add the 
probeID of the bead type to a list which we return.

The standard deviation that we calculate is based on the log2 bead intensities.
}

\value{

numeric vector giving the probeIDs which were calculated to have standard deviation greater than 'limit'.

}

\author{Mark Dunning}

\examples{

data(BLData)

##Default is to find bead types with standard deviation of intensity > 1

findHighestSD(BLData, 1)

##limit can be changed eg. to give ProbeIDs with standard deviation > 0.7

findHighestSD(BLData, 1, limit=0.7)

##Any ProbeIDs can be examined in more detail

plotBeadIntensities(BLData, 678,1)

probeDiagnostics(BLData, 678, 1)

}

\keyword{methods}
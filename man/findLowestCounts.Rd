\name{findLowestCounts}
\alias{findLowestCounts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find Lowest bead counts}
\description{
Find bead types which are represented less than 24 times a given array
}
\usage{
findLowestCounts(BLData, array, limit = 24)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{BeadLevelList object containing bead level data}
  \item{array}{numeric value for the array we are interested in}
  \item{limit}{numeric value defining when bead types should be returned by the function}
}
\details{
It is known that the 
numbers of each bead type on an array should be Poisson and on average around 30 of each bead type should be seen. This high amount of replication is 
a major advantage of the BeadArray technology.

The function simply loops through all the different ProbeIDs found on the array and counts how many of the bead type there are. If the number of beads found 
is less than 'limit' then the bead type is added to a list which is returned by the function.

The default value for 'limit' is 24 which is the 5th percentile of the appropriate Poisson distribution. 

WARNING: THIS FUNCTION WILL NOT BE EFFICIENT FOR WHOLE GENOME DATA
}
\value{
  numeric vector giving the ProbeIDs of the bead types which occur less than 'limit' times on the array.
}

\author{Mark Dunning}

\examples{

data(BLData)

##Default is to find all bead types represented less than 24 times

findLowestCounts(BLData, 1)

##Can change limit parameter to give all bead types represented less than 20 times

findLowestCounts(BLData, 1, limit=20)

##ProbeIDs returned by function can be investigated in more detail

plotBeadLocations(BLData, ProbeID=87, array=1)

}
\keyword{methods}




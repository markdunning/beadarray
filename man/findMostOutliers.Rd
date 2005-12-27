\name{findMostOutliers}
\alias{findMostOutliers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Beads with most outliers}
\description{
Function to find out which bead types on a particular array have the most outliers}


\usage{
findMostOutliers(BLData, array, limit = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{an BeadLevelList object containing bead level data}
  \item{array}{numeric value giving the number of the array we want to find outliers for} 
  \item{limit}{numeric value giving the limit for the number of outliers (default is 5)}
}
\details{
Takes each bead type found on the array in turn and calculates the number of outliers found for the bead type. If more than 'limit' outliers are found for a bead type we add the probeID
to a list. Preliminary investigation suggests that bead types rarely have more than 5 outliers.

This function uses the default method of finding outliers using 3 median absolute deviations from the bead mean with raw (un-logged data).

}
\value{
vector of numeric values giving the probeIDs of bead types which have more than 'limit' number of outliers.
}

\author{Mark Dunning}
\seealso{ \code{\link{findOutliers}} }

\examples{

data(BLData)

##Default is to find bead types with more than 5 outliers

findMostOutliers(BLData, 1)

##Can change limit parameter to find all bead types with more than 3 outliers (say)

findMostOutliers(BLData, 1, limit=3)

##ProbeIDs returned by function can be investigated in more detail

plotBeadIntensities(BLData, 807, 1)


}

\keyword{methods}
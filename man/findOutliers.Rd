\name{findOutliers}
\alias{findOutliers}
\alias{getProbeIndices}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Find Outliers }
\description{
 Function finds all beads which are outliers for a given bead type
}
\usage{
findOutliers(BLData, probe, array, log=FALSE, n=3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{BeadLevelList object containing bead level data}
  \item{probe}{numeric value for the probeID of the bead type we want to find outliers for}
  \item{array}{numeric value for the array we want to find outliers on}
  \item{log}{if TRUE the log2 intensities of each bead will be used to find outliers. Otherwise un-logged data is used}
  \item{n}{numeric value defining a cut-off of 'n' median absolute deviations from the mean to determine outliers. The default value is 3}
}
\details{

The intensities of each bead with probeID 'probe' on the specified array are found and if the 'log' parameter is set to true we do a log2 transformation
on the set of values.

The mean and median absolute deviation of all the bead intensities are then calculated. Outliers are beads which have intensity more than 
'n' median absolute deviations from the mean.


The method used by Illumina is to use un-logged intensities with a value 'n' = 3.

Any beads which have intensity NA are also counted as outliers

}

\value{
  
  List of beadIDs (1-49777) dividing the beads of this bead type into two categories.

  \item{valid}{valid beads}
  \item{outliers}{beads which are calculated as outliers}
}

\author{Mark Dunning}

\seealso{ \code{\link{findAllOutliers}} }

\examples{

data(BLData)

findOutliers(BLData, 2, 1)

findOutliers(BLData, 2, 1, log=TRUE)


##Methods may not always agree

findOutliers(BLData, 6140, 1)

findOutliers(BLData, 6140, 1, log=TRUE)

##Increasing value of n selects less outliers

findOutliers(BLData, 6140, 1, n=4)

findOutliers(BLData, 6140, 1, log=TRUE, n=4)

}

\keyword{methods}
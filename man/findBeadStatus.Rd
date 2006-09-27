\name{findBeadStatus}
\alias{findBeadStatus}
\alias{getProbeIndices}
\alias{getProbeIndicesC}
\alias{getProbeIntensities,BeadLevelList-method}
\alias{getProbeIntensities}
\alias{findAllOutliers, BeadLevelList-method}
\alias{createBeadSummaryData, BeadLevelList-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Find Outliers }
\description{
 Function finds all beads which are outliers for a given bead type
}
\usage{
findBeadStatus(BLData, probes, array, log=FALSE, n=3, outputValid = FALSE, intProbeID = NULL, ignoreList=NULL, probeIndex=NULL, startSearch=1)
getProbeIntensities(BLData, ProbeIDs, array,log=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{BeadLevelList object containing bead level data}
  \item{probes}{numeric vector for the ProbeIDs of the bead type we want to find outliers for}
  \item{array}{numeric value for the array we want to find outliers on}
  \item{log}{if TRUE the log2 intensities of each bead will be used to find outliers. Otherwise un-logged data is used}
  \item{n}{numeric value defining a cut-off of 'n' median absolute deviations from the mean to determine outliers. The default value is 3}
  \item{outputValid}{if TRUE the IDs of beads which are not outliers will be output}
  \item{intProbeID}{BLData\$ProbeID coerced to vector of integers.  Never change this, for internal use only}
\item{ignoreList}{list of ProbeIDs to be omitted from the averaging procedure. These could be Illumina internal controls which are replicated many thousands of times on arrays}
\item{probeIndex}{parameter for internal use only}
\item{startSearch}{integer specify where to start searching for a particular ProbeID}
\item{ProbeIDs}{numeric vector for the ProbeIDs of the bead type we want to find outliers for}

}
\details{

The intensities of each bead with ProbeID 'probe' on the specified array are found and if the 'log' parameter is set to true we do a log2 transformation
on the set of values.

The mean and median absolute deviation of all the bead intensities are then calculated. Outliers are beads which have intensity more than 
'n' median absolute deviations from the mean.


The method used by Illumina is to use un-logged intensities with a value 'n' = 3.

Any beads which have intensity NA are also counted as outliers

The function returns only the outliers for a bead type unless the
outputValid parameter is specified.

}

\value{
  
  List of beadIDs (1-49777) dividing the beads of this bead type into two categories.

  \item{valid}{valid beads}
  \item{outliers}{beads which are calculated as outliers}
}

\author{Mike Smith and Mark Dunning}

\seealso{ \code{findAllOutliers} }

\examples{

data(BLData)

findBeadStatus(BLData, 70450, 1, outputValid=TRUE)

findBeadStatus(BLData, 70450, 1, log=TRUE, outputValid=TRUE)

getProbeIntensities(BLData, 70450, 1)



}

\keyword{methods}
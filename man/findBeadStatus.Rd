\name{findBeadStatus}
\alias{findBeadStatus}
\alias{getProbeIntensities}
%\alias{getProbeIndices}
%\alias{getProbeIndicesC}
\alias{getProbeIntensities,BeadLevelList-method}
%\alias{findAllOutliers, BeadLevelList-method}
%\alias{createBeadSummaryData, BeadLevelList-method}

\title{ Find Outliers}
\description{
 Function finds all beads which are outliers for a given bead type
}
\usage{
findBeadStatus(BLData, probes, array = 1, log = FALSE, what = "G", n = 3,
           outputValid = FALSE, intProbeID = NULL, ignoreList = NULL,
           probeIndex = NULL, startSearch = 1)
getProbeIntensities(BLData, ProbeIDs, array = 1, log = TRUE, what = "G")
}
\arguments{
  \item{BLData}{\code{BeadLevelList}}
  \item{probes}{numeric vector for the ProbeIDs of the bead type we
    want to find outliers for}
  \item{array}{integer specifying which strip/array to use}
  \item{log}{if TRUE the intensities will be calculated on the 
    log2 scale. Otherwise un-logged data is used}
  \item{what}{character string specifying which intensities to use.
    Possibilities are \code{"G"}, \code{"Gb"} for single channel data
    and \code{"G"}, \code{"Gb"}, \code{"R"} and \code{"Rb"} for
    two-colour data}
  \item{n}{numeric value defining a cut-off for the number of 
    median absolute deviations (MADs) from the median to use for
    determining outliers.  The default value is 3.}
  \item{outputValid}{if TRUE the IDs of beads which are not outliers
    will be output}
  \item{intProbeID}{BLData\$ProbeID coerced to vector of integers.
    Never change this, for internal use only}
  \item{ignoreList}{list of ProbeIDs to be omitted from the averaging
    procedure. These could be Illumina internal controls which are
    replicated many thousands of times on arrays}
  \item{probeIndex}{parameter for internal use only}
  \item{startSearch}{integer specify where to start searching for a
    particular ProbeID}
  \item{ProbeIDs}{numeric vector for the ProbeIDs of the bead type
    we want to find outliers for}
}
\details{
  The intensities of each bead with ProbeID 'probe' on the
  specified array are found and if the 'log' parameter is set
  to TRUE we do a log2 transformation of these values.

  The median and MAD for the bead intensities are then calculated. 
  Outliers are beads which have intensity more
  than 'n' MADs from the median.

  The method used by Illumina is to use un-logged intensities with
  'n' = 3.

  Any beads which have intensity NA are also counted as outliers.

  The function returns only the outliers for a bead type unless the
  outputValid parameter is specified.
}

\value{
  List of beadIDs dividing the beads of this bead type into two
  categories.
  \item{valid}{valid beads}
  \item{outliers}{beads which are calculated as outliers}
}

\author{Mike Smith and Mark Dunning}

\seealso{ \code{findAllOutliers} }

\examples{
data(BLData)
findBeadStatus(BLData, 2, 1, outputValid=TRUE)
findBeadStatus(BLData, 2, 1, log=TRUE, outputValid=TRUE)
findBeadStatus(BLData, 23, 1, outputValid=TRUE)
findBeadStatus(BLData, 23, 1, log=TRUE, outputValid=TRUE)
}
\keyword{methods}

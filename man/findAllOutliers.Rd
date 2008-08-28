\name{findAllOutliers}
\alias{findAllOutliers}
%\alias{findAllOutliers-methods}
\title{Find outliers on a given strip/array}
\description{
  Function to find all beads which are outliers for their
  particular bead type on a given strip/array using Illumina's 
  standard outlier detection method.
}
\usage{
findAllOutliers(BLData, array, log=FALSE, n=3, what="G", usewts=FALSE)
}
\arguments{
  \item{BLData}{\code{BeadLevelList}}
  \item{array}{integer specifying which strip/array we
    want to find outliers on}
  \item{log}{if TRUE the intensities will be calculated on the log2 scale. Otherwise un-logged data is used.}
  \item{n}{numeric value defining a cut-off for the number of 
    median absolute deviations (MADs) from the median to use for 
    determining outliers.  The default value is 3.}
  \item{what}{character string specifying which intensities to use.
    See \code{\link{getArrayData}} for a list of possibilities.}
  \item{usewts}{if TRUE, then beads with weights below 1 will be discarded prior to analysis.}
}
\details{
  We find the outliers for each bead type on the array in
  turn using the \code{\link{findBeadStatus}} function and
  store the indices of the outliers found. By default, outliers
  for a particular bead type are determined using a 3
  MAD cut-off from the median.
}
\value{
  numeric vector giving the row indices of \code{BLData} (in the
  range 1 to total number of beads on the array) of all beads 
  that are outliers for their bead type.
}
\author{Mark Dunning}
\seealso{\code{\link{findBeadStatus}}}
\examples{
data(BLData)
# how many outliers are there on the original scale?
length(findAllOutliers(BLData, 1))
# how many outliers are there on the log2-scale?
length(findAllOutliers(BLData, 1, log=TRUE)) #
# how many outliers are there using a 4 MAD 
# cut-off from the median?
length(findAllOutliers(BLData, 1, n=4))
}
\keyword{methods}

\name{findAllOutliers}
\alias{findAllOutliers}
\alias{findAllOutliers-methods}
\title{Find outliers on array}
\description{
  Function to find all beads which are outliers for their
  particular bead type on a given array.
}

\usage{
findAllOutliers(BLData, array, log=FALSE, n=3, which="G")
}
\arguments{
  \item{BLData}{a BeadLevelList object}
  \item{array}{numeric value specifying which array we
    want to find outliers on}
  \item{log}{if TRUE the log2 intensities will be
    used to calculate outliers}
  \item{n}{Specify a cut-off for outliers as n median
    absolute deviations from the median. The default value is 3}
  \item{which}{character string specifying which intensities to plot.
    Possibilities are \code{"G"}, \code{"Gb"} for single channel data
    and \code{"G"}, \code{"Gb"}, \code{"R"} and \code{"Rb"} for
    two-colour data}
}
\details{
  We find the outliers for each bead type on the array in
  turn using the \code{\link{findBeadStatus}} function and
  keep a list of the outliers found. By default, outliers
  for a particular bead type are determined using a 3
  MAD cut-off from the mean.
}
\value{
  numeric vector giving the row indices of BLData (in the
  range 1 to total number of beads on the array) of all
  the beads that are outliers for their bead type.
}
\author{Mark Dunning}
\seealso{\code{\link{findBeadStatus}}}
\examples{
data(BLData)

length(findAllOutliers(BLData, 1))

length(findAllOutliers(BLData, 1, log=TRUE))

length(findAllOutliers(BLData, 1, n=5))
}
\keyword{methods}
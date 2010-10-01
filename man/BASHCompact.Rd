\name{BASHCompact}
\alias{BASHCompact}
\title{BASH - Compact Defect Analysis}
\description{
Creates a list of probes marked as being in compact defects.
}
\usage{
BASHCompact(BLData, array, neighbours = NULL, useLocs = TRUE, transFun = logGreenChannelTransform, maxiter = 10, cutoff = 8, cinvasions = 10, ...)
}

\arguments{
  \item{BLData}{\code{BeadLevelList}}
  \item{array}{integer specifying which strip/array to plot}
  \item{neighbours}{A Neighbours matrix. Optional - if left NULL, it will be computed.}
  \item{useLocs}{logical value, specifiying whether the .locs file (if present) should be used to determine neighbours.}
  \item{transFun}{function to use to transform data prior to running BASH.}
  \item{maxiter}{Integer - Maximum number of iterations.}
  \item{cutoff}{Integer - Size a cluster must be to be labelled a compact defect.}
  \item{cinvasions}{Integer - Number of invasions used when closing the image.}
  \item{...}{Additional arguments to be passed to \code{findAllOutliers}.}
}

\details{
	\code{BASHCompact} finds "compact defects" on an array. A compact defect is defined as a large connected cluster of outliers.

	This function first finds the outliers on an array. This is done via the function \code{findAllOutliers}.

	Next, using the Neighbours matrix and a Flood Fill algorithm, it determines which beads are in large connected clusters of outliers (of size larger than \code{cutoff}). These beads are then temporarily removed and the process repeated with the remaining beads. The repetition continues until either no large clusters of outliers remain, or until we have repeated the process \code{maxiter} times (and in this case, a warning will be given). In this way, we obtain a list of defective probes.

	Finally, we "close" the image, to fill in small gaps in the defect image. This consists of a "dilation" and an "erosion". In the dilation, we expand the defect image, by adding beads adjacent to defective beads into the defect image. This is repeated \code{cinvasions} times. In the erosion, we contract the defect image, by removing beads adjacent to non-defective beads from the defect image. (Erosion of the defect image is equivalent to a dilation of the non-defective image.)
}

\value{
A vector consisting of the BeadIDs of beads labelled as compact defects.
}

\author{Jonathan Cairns}

\seealso{\code{\link{BASH}}, \code{\link{generateNeighbours}} }

\examples{
data(BLData)
o <- BASHCompact(BLData, 1)
o <- BASHCompact(BLData, 1, cinvasions = 10) ##increased no of closure invasions
o <- BASHCompact(BLData, 1, cutoff = 12) ##only larger defects will be found with this setting
}
\references{
	Mayte Suarez-Farinas, Maurizio Pellegrino, Knut M. Wittkwosky and Marcelo O. Magnasco (2007). Harshlight: A "corrective make-up" program for microarray chips. R package version 1.8.0. http://asterion.rockefeller.edu/Harshlight/

}
\keyword{misc}

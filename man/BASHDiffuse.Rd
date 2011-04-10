\name{BASHDiffuse}
\alias{BASHDiffuse}
\title{BASH - Diffuse Defect Analysis}
\description{
Creates a list of probes marked as being in diffuse defects.
}
\usage{
BASHDiffuse(BLData, array, transFun = logGreenChannelTransform, neighbours = NULL, useLocs = TRUE, E = NULL, n = 3, compact = NULL, sig = 0.0001, invasions = 10, cutoff = 8, cinvasions = 10, twotail = FALSE)
}

\arguments{
  \item{BLData}{\code{BeadLevelList}}
  \item{array}{integer specifying which section/array to analyse}
  \item{transFun}{Function to use to transform data prior to running BASHDiffuse.}
  \item{neighbours}{A Neighbours matrix. Optional - if left NULL, it will be computed, using default \code{\link{generateNeighbours}} settings.}
  \item{useLocs}{Logical value, specifiying whether the .locs file (if present) should be used to determine neighbours.}
  \item{E}{Numerical vector - The error image to use. Optional - if left blank, it will be computed, using \code{generateE} using \code{bgfilter = "median"}.}
  \item{n}{Specify a cut-off for outliers as n median absolute deviations (MADs) from the median. The default value is 3}
  \item{compact}{Vector - Optional. BeadIDs of beads in compact defects to remove from the analysis.}
  \item{sig}{Numerical - Significance level of binomial test.}
  \item{invasions}{Integer - Number of invasions to use to find the kernel (see below).}
  \item{cutoff}{Integer - Size a cluster must be to be labelled a diffuse defect.}
  \item{cinvasions}{Integer - Number of invasions used when closing the image.}
  \item{twotail}{Logical - If TRUE, then we analyse positive and negative outliers separately, and then combine the diffuse defect images at the end.}
}

\details{
	\code{BASHDiffuse} finds "diffuse defects" on an array. A diffuse defect is defined as a region containing an unusually large number of (not necessarily connected) outliers.

	Firstly, we consider the error image \code{E}, and find outlier beads on this image. Outliers for a particular bead type are determined using a 3 MAD cut-off from the median.

	We now consider an area around each bead (known as the "kernel"). The kernel is found by an invasion process using the neighbours matrix - we choose the beads which can be reached from the central bead in \code{cinvasions} steps.

	We count how many beads are in the kernel, and how many of these are marked as outliers. Using a binomial test, we work out if there are significantly more outliers in the kernel than would be expected if the outliers were equally distributed over the entire array. If so, then the central bead is marked as a diffuse defect.

	Lastly, we run a clustering algorithm and a closing algorithm similar to those in \code{\link{BASHCompact}}.
}

\value{
A vector consisting of the BeadIDs of beads considered diffuse defects.
}

\author{Jonathan Cairns}

\seealso{\code{\link{BASH}}, \code{\link{generateNeighbours}}, }

\examples{
data(BLData)
o <- BASHDiffuse(BLData, 1)
o <- BASHDiffuse(BLData, 1, sig = 0.00001) ##stricter significance value, perhaps more useful on a BeadChip.
o <- BASHDiffuse(BLData, 1, cutoff = 12) ##only larger defects will be found with this setting
}
\references{
	Mayte Suarez-Farinas, Maurizio Pellegrino, Knut M. Wittkwosky and Marcelo O. Magnasco (2007). Harshlight: A "corrective make-up" program for microarray chips. R package version 1.8.0. http://asterion.rockefeller.edu/Harshlight/

}
\keyword{misc}

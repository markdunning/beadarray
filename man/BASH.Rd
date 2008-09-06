\name{BASH}
\alias{BASH}
\title{BASH - BeadArray Artefact Selection, using Harshlight}
\description{
BASH is an automatic detector of physical defects on an array. It is designed to detect three types of defect - COMPACT, DIFFUSE and EXTENDED.
}
\usage{
BASH(BLData, array = NULL, compact = TRUE, diffuse = TRUE, extended = TRUE, log = TRUE, cinvasions = 10, dinvasions = 15, einvasions = 20, bgcorr = "median", maxiter = 10, compcutoff = 8, compdiscard = TRUE, diffcutoff = 10, diffsig = 0.0001, diffn = 3, difftwotail = FALSE)
}
\arguments{
  \item{BLData}{\code{BeadLevelList}}
  \item{array}{integer specifying which strip/array to plot. Alternatively you can supply a vector of strip/array IDs, and BASH will analyse each in turn.}
  \item{compact}{Logical - Perform compact analysis?}
  \item{diffuse}{Logical - Perform diffuse analysis?}
  \item{extended}{Logical - Perform extended analysis?}
  \item{log}{Logical - Perform analyses on the log scale? (recommended)}
  \item{cinvasions}{Integer - number of invasions used whenever closing the image - see \code{\link{BASHCompact}}}
  \item{dinvasions}{Integer - number of invasions used in diffuse analysis, to find the kernel - see \code{\link{BASHDiffuse}}}
  \item{einvasions}{Integer - number of invasions used when filtering the error image - see \code{\link{BGFilter}}.}
  \item{bgcorr}{One of "none", "median", "medianMAD" - Used in diffuse analysis, this determines how we attempt to compensate for the background varying across an array. For example, on a SAM array this should be left at "median", or maybe even switched to "none", but if analysing a large beadchip then you might consider setting this to "medianMAD". (this code is passed to the \code{method} argument of \code{\link{BGFilter}}).}
  \item{maxiter}{Integer - Used in compact analysis - the max number of iterations allowed. (Exceeding this results in a warning.)}
  \item{compcutoff}{Integer - the threshold used to determine whether a group of outliers is in a compact defect. In other words, if a group of at least this many connected outliers is found, then it is labelled as a compact defect.}
  \item{compdiscard}{Logical - should we discard compact defect beads before doing the diffuse analyis?}
  \item{diffcutoff}{Integer - this is the threshold used to determine the minimum size that clusters of diffuse defects must be.}
  \item{diffsig}{Probability - The significance level of the binomial test performed in the diffuse analysis.}
  \item{diffn}{Numerical - when finding outliers on the diffuse error image, how many MADs away from the median an intensity must be for it to be labelled an outlier.}
  \item{difftwotail}{Logical - If TRUE, then in the diffuse analysis, we consider the high outlier and low outlier images seperately.}
}
\details{
	The \code{BASH} pipeline function performs three types of defect analysis on an image.

	The first, COMPACT DEFECTS, finds large clusters of outliers, as per \code{BASHCompact}. The outliers are found using \code{findAllOutliers()}. We then find which outliers are clustered together. This process is iterative - having found a compact defect, we remove it, and then see if any more defects are found.

	The second, DIFFUSE DEFECTS, finds areas which are densely populated with outliers (which are not necessarily connected), as per \code{BASHDiffuse}. To make this type of defect more obvious, we first generate an ERROR IMAGE, and then find outliers based on this image. (The error image is calculated by using \code{method = "median"} and \code{bgfilter = "medianMAD"} in \code{generateE}, unless \code{ebgcorr = FALSE} in which case we use \code{bgfilter = "median"}.) Now we consider a neighbourhood around each bead and count the number of outlier beads in this region. Using a binomial test we determine whether this is more that we would expect if the outliers were evenly spread over the entire array. If so, we mark it as a diffuse defect. (A clustering algorithm similar to the compact defect analysis is run to reduce false positives.)

	After each of these two analyses, we "close" the image, filling in gaps.

	The third, EXTENDED DEFECTS, returns a score estimating how much the background is changing across an array, as per \code{BASHExtended}. To estimate the background intensity, we generate an error image using the median filter (i.e. \code{generateE} with \code{method = "median"} and \code{bgfilter = "median"}). We divide the variance of this by the variance of an error image without using the median filter, to obtain our extended score.

	It should be noted that to avoid repeated computation of distance, a "neighbours" matrix is used in the analysis. This matrix describes which beads are close to other beads. If a large number of beads are missing (for example, if beads with ProbeID = 0 were discarded) then this algorithm may be affected.

	For more detailed descriptions of the algorithms, read the help files of the respective functions listed in "see also".
}

\value{
	The output has three attributes:
	
	$wts: A list, where the ith object in the list corresponds to the weights for array i.

	$ext: A vector of extended scores (null if the extended analysis was disabled)

	$call: The function you used to call BASH.
}

\author{Jonathan Cairns}

\seealso{\code{\link{BASHCompact}}, \code{\link{BASHDiffuse}}, \code{\link{BASHExtended}}, \code{\link{generateE}}, \code{\link{generateNeighbours}}, }

\examples{
	data(BLData)
	output <- BASH(BLData)
	boxplot(output$ext) #view spread of extended scores
	for(i in 1:4)
	{
		BLData <- setWeights(BLData, output$wts[[i]], i) #apply BASH weights to BLData
	}

	#diffuse test is stricter
	output <- BASH(BLData, diffsig = 0.00001)

	#more outliers on the error image are used in the diffuse analysis
	output <- BASH(BLData, diffn = 2)

	#only perform compact & diffuse analyses (we will only get weights)
	output <- BASH(BLData, extended = FALSE)

	#attempt to correct for background.
	output <- BASH(BLData, ebgcorr = TRUE)
}
\references{
	Mayte Suarez-Farinas, Maurizio Pellegrino, Knut M. Wittkwosky and Marcelo O. Magnasco (2007). Harshlight: A "corrective make-up" program for microarray chips. R package version 1.8.0. http://asterion.rockefeller.edu/Harshlight/

	#FIXME
}
\keyword{misc}

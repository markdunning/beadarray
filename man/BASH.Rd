\name{BASH}
\alias{BASH}

\title{BASH - BeadArray Subversion of Harshlight}
\description{
BASH is an automatic detector of physical defects on an array. It is designed to detect three types of defect - COMPACT, DIFFUSE and EXTENDED.
}
\usage{
BASH(BLData, array, neighbours=NULL, transFun = logGreenChannelTransform,
    outlierFun = illuminaOutlierMethod, compn=3, wtsname=NULL, compact = TRUE, 
    diffuse = TRUE, extended = TRUE, cinvasions = 10, dinvasions = 15, 
    einvasions = 20, bgcorr = "median", maxiter = 10, compcutoff = 8, 
    compdiscard = TRUE, diffcutoff = 10, diffsig = 0.0001, diffn = 3, 
    difftwotail = FALSE, useLocs = TRUE, ...)
}
\arguments{
  \item{BLData}{\code{BeadLevelList}}
  \item{array}{integer specifying which section/array to plot.}
  \item{neighbours}{the user may specify the neighbours matrix, rather than have BASH calculate it. Time can be saved if using BASH and \code{\link{HULK}}, by calculating the neighbours matrix once and passing it to the two functions.}
  \item{transFun}{function to use to transform data prior to running BASH}
  \item{outlierFun}{the choice of outlier calling function to use.}
   \item{compn}{Numerical - when finding outliers in the compact analysis, how many MADs away from the median (for example) an intensity must be for it to be labelled an outlier.}
  \item{wtsname}{name under which bead weights are stored in the BLData object. It is only necessary to specify this if a) weights have already been set, and b) you wish BASH to observe them.}
  \item{compact}{Logical - Perform compact analysis?}
  \item{diffuse}{Logical - Perform diffuse analysis?}
  \item{extended}{Logical - Perform extended analysis?}
  \item{cinvasions}{Integer - number of invasions used whenever closing the image - see \code{\link{BASHCompact}}}
  \item{dinvasions}{Integer - number of invasions used in diffuse analysis, to find the kernel - see \code{\link{BASHDiffuse}}}
  \item{einvasions}{Integer - number of invasions used when filtering the error image - see \code{BGFilter}.}
  \item{bgcorr}{One of "none", "median", "medianMAD" - Used in diffuse analysis, this determines how we attempt to compensate for the background varying across an array. For example, on a SAM array this should be left at "median", or maybe even switched to "none", but if analysing a large beadchip then you might consider setting this to "medianMAD". (this code is passed to the \code{method} argument of \code{BGFilter}). Note that "none" may be the correct setting if \code{HULK} has already been applied.}
  \item{maxiter}{Integer - Used in compact analysis - the max number of iterations allowed. (Exceeding this results in a warning.)}
  \item{compcutoff}{Integer - the threshold used to determine whether a group of outliers is in a compact defect. In other words, if a group of at least this many connected outliers is found, then it is labelled as a compact defect.}
  \item{compdiscard}{Logical - should we discard compact defect beads before doing the diffuse analyis?}
  \item{diffcutoff}{Integer - this is the threshold used to determine the minimum size that clusters of diffuse defects must be.}
  \item{diffsig}{Probability - The significance level of the binomial test performed in the diffuse analysis.}
  \item{diffn}{Numerical - when finding outliers on the diffuse error image, how many MADs away from the median an intensity must be for it to be labelled an outlier.}
  \item{difftwotail}{Logical - If TRUE, then in the diffuse analysis, we consider the high outlier and low outlier images seperately.}
  \item{useLocs}{Logical - If TRUE then a .locs file corresponding to the array is sought and, if found, used to identify the neighbouring beads.  If FALSE the neighbours are infered algorithmically.  See \code{\link{generateNeighbours}} for more details.}
  \item{...}{Logical - Perform compact analysis?}
}
\details{
	The \code{BASH} pipeline function performs three types of defect analysis on an image.

	The first, COMPACT DEFECTS, finds large clusters of outliers, as per \code{BASHCompact}. The outliers are found using \code{findAllOutliers()}. We then find which outliers are clustered together. This process is iterative - having found a compact defect, we remove it, and then see if any more defects are found.

	The second, DIFFUSE DEFECTS, finds areas which are densely populated with outliers (which are not necessarily connected), as per \code{BASHDiffuse}. To make this type of defect more obvious, we first generate an ERROR IMAGE, and then find outliers based on this image. (The error image is calculated by using \code{method = "median"} and \code{bgfilter = "medianMAD"} in \code{generateE}, unless \code{ebgcorr = FALSE} in which case we use \code{bgfilter = "median"}.) Now we consider a neighbourhood around each bead and count the number of outlier beads in this region. Using a binomial test we determine whether this is more that we would expect if the outliers were evenly spread over the entire array. If so, we mark it as a diffuse defect. (A clustering algorithm similar to the compact defect analysis is run to reduce false positives.)

	After each of these two analyses, we "close" the image, filling in gaps.

	The third, EXTENDED DEFECTS, returns a score estimating how much the background is changing across an array, as per \code{BASHExtended}. To estimate the background intensity, we generate an error image using the median filter (i.e. \code{generateE} with \code{method = "median"} and \code{bgfilter = "median"}). We divide the variance of this by the variance of an error image without using the median filter, to obtain our extended score.

	It should be noted that to avoid repeated computation of distance, a "neighbours" matrix is used in the analysis. This matrix describes which beads are close to other beads. If a large number of beads are missing (for example, if beads with ProbeID = 0 were discarded) then this algorithm may be affected.

	For more detailed descriptions of the algorithms, read the help files of the respective functions listed in "see also".

	BASH is currently quite a slow, memory-intensive function. It will only run on a single array at a time, and for analysis of multiple arrays, we recommend parallelising the command. An example is shown using the base parallel package.
      
}

\value{
	The output is a list with four attributes:
	
	wts: A vector of weights for the matrix.

	ext: A vector of extended scores (null if the extended analysis was disabled).

  QC: A summary of the extended score and the number of beads masked.

	call: The function you used to call BASH.
}

\author{Jonathan Cairns}

\seealso{\code{\link{BASHCompact}}, \code{\link{BASHDiffuse}}, \code{\link{BASHExtended}}, \code{\link{generateNeighbours}}, \code{\link{HULK}} }

\examples{

\dontrun{

if(require(beadarrayExampleData)){

	data(exampleBLData)
	output <- BASH(exampleBLData,array=1,useLocs=FALSE)
        exampleBLData <- setWeights(exampleBLData, output$wts, array=1) #apply BASH weights to exampleBLData
	
	###BASH only accepts one array at a time, but it can be made to run in a parallel fashion
	library(parallel)

	output <- mclapply(c(1,2), function(x) BASH(exampleBLData, array=x, useLocs=FALSE))

	for(i in 1:2){
	  exampleBLData <- setWeights(exampleBLData, output[[i]]$wts, array=i) 
	}
      

	#diffuse test is stricter
	output <- BASH(exampleBLData, diffsig = 0.00001,array=1, useLocs=FALSE)

	#more outliers on the error image are used in the diffuse analysis
	output <- BASH(exampleBLData, diffn = 2,array=1, useLocs=FALSE)

	#only perform compact & diffuse analyses (we will only get weights)
	output <- BASH(exampleBLData, extended = FALSE,array=1, useLocs=FALSE)

	#attempt to correct for background.
	output <- BASH(exampleBLData, bgcorr = "median",array=1, useLocs=FALSE)
}

else{
  
  stop("You will need the beadarrayExampleData package to run this example")
}



}


}

\references{
	J. M. Cairns, M. J. Dunning, M. E. Ritchie, R. Russell, and A. G. Lynch (2008). BASH: a tool for managing BeadArray spatial artefacts. Bioinformatics 15; 24(24)

}
\keyword{misc}

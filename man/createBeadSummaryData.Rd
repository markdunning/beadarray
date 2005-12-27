\name{createBeadSummaryData}
\alias{createBeadSummaryData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Produce bead averages}
\description{
Produce bead averages for each bead type used in an experiment on a specified set of arrays.
}
\usage{
createBeadSummaryData(BLData, log=FALSE, n=3, arrays=seq(1:length(BLData$R[1,])))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{an BeadLevelList object containing bead level data}
  \item{log}{if TRUE then use the log2 intensities of each bead when removing outliers}
  \item{n}{use n median absolute deviations from the mean as a cut-off for outliers} 	
  \item{arrays}{numeric vector giving the indices of the arrays we want to average. Default is all arrays in the BeadLevelList}
}
\details{

'arrays' specifies which arrays in the BeadLevelList object we want to average. The default is to average all arrays in the BeadLevelList object but we might only want to average the 
half of the arrays to save time. Objects which are created seperately by 'createBeadSummaryData' may be joined by using the \code{\link{cbind}} function. 

To summarise each bead type we first have to remove
the outliers for the bead type on the array. The default method is the method used by Illumina which uses the raw (un-logged) intensities and categorises outliers as being those beads 
which have intensity more than 3 median absolute deviations from the bead mean intensity. By setting the log parameter we can use the log2 intensities instead
of un-logged and by setting the 'n' parameter we choose beads with intensity 'n' mads from the mean.


With outliers removed, we simply take an average of the red and green intensities of 
the remaining beads and create \$R and \$G matrices. We also store the number of outliers that were detected for each bead type and the standard deviation of the 
remaining beads in the \code{$nooutliers} and \code{$beadstdev} columns respectively.

As the function runs it prints out the index of the array that is currently being processed.

}
\value{
  An BeadSummaryList (bead summary) object in which all components (R, Rb, G, ....) are matrices with number of rows equal to the number of bead types for the experiment and number of columns equal to the length of the 'arrays' object.
}

\author{Mark Dunning}



\seealso{
\code{\link{findOutliers}}	

}

\examples{

##produce bead summaries for the first array in our BLData object

#data(BLData)
#BSData = createBeadSummaryData(BLData)
#names(BSData)

}

\keyword{methods}
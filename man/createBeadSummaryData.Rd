\name{createBeadSummaryData}
\alias{createBeadSummaryData,BeadLevelList-method}
\alias{createBeadSummaryData}
\title{Produce bead averages}
\description{
  Produce bead averages for each bead type used in an experiment
  on a specified set of arrays.
}
\usage{
createBeadSummaryData(BLData, log=FALSE, n=3,
               imagesPerArray = 1, what="G",  probes = NULL)
}
\arguments{
  \item{BLData}{an BeadLevelList object}
  \item{log}{if TRUE then summarise the log2 intensities of each bead}
  \item{n}{use n median absolute deviations from the median as a
    cut-off for outliers}
  \item{imagesPerArray}{Specifies how many images (strips) there are
    per array.  Normally 1 for a SAM and 1 or 2 for a BeadChip.
    The images (strips) from the same array will be combined so that
    each column in the output represents a sample
  }
  \item{what}{character string specifying which intensities/values
    to summarise. Choices are \code{"G"}, \code{"R"}, \code{"RG"},
    \code{"M"} or  \code{"A"}}
  \item{probes}{Specify particular probes to summarise. If left NULL
    then all the probes on the first array are used.}
}
\details{
  Objects which are created seperately by 'createBeadSummaryData' may
  be joined by using the \code{\link{cbind}} function.

  To summarise each bead type we first have to remove the outliers for
  the bead type on the array. The default method is the method used by
  Illumina which uses the raw (un-logged) intensities and categorises
  outliers as being those beads which have intensity more than 3 
  median absolute deviations from the bead median intensity. By setting
  the log argument to \code{TRUE} we calculate outliers and summary 
  values on the log2-scale.  By changing the 'n' argument we choose
  beads with intensity 'n' mads from the median.

  With outliers removed, we calculate the average (mean) intensities
  of the remaining beads as well as the standard error and number of beads.
}
\value{
  An ExpressionSetIllumina (bead summary) object (or SnpSetIllumina
  object for two-colour data, when \code{what="RG"}) in which all
  components are matrices with number of rows equal to the number of bead
  types for the experiment and number of columns equal to the number
  of arrays.
}
\author{Mark Dunning and Mike Smith}
\seealso{
\code{\link{findBeadStatus}}	
}
\examples{
#produce bead summaries for the first array in our BLData object
data(BLData)
BSData = createBeadSummaryData(BLData, log=TRUE)
dim(BSData)
}
\keyword{methods}

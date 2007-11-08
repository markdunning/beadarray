\name{createBeadSummaryData}
\alias{createBeadSummaryData,BeadLevelList-method}
\alias{createBeadSummaryData}
\title{Produce bead averages}
\description{
  Produce bead averages for each bead type used in an experiment
  on a specified set of arrays.
}
\usage{
createBeadSummaryData(BLData, log=FALSE, imagesPerArray = 1,
               what="G",  probes = NULL, arrays=NULL, 
               method="illumina", n=3, trim=0.05)
}
\arguments{
  \item{BLData}{an BeadLevelList object}
  \item{log}{if TRUE then summarise the log2 intensities of each bead}
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
  \item{arrays}{Specify particular arrays to summarise. If left NULL
    then all arrays are summarised.}
  \item{method}{chracter string specifying the summarisation method to
    use. Options are \code{"illumina"}, \code{"mean"}, \code{"median"},
    \code{"trim"} and \code{"winsorise"}.}
  \item{n}{use n median absolute deviations from the median as a
    cut-off for outliers}
  \item{trim}{fraction of intensities to remove from the bead summary
    calculations when \code{method="trim"}, or the fraction of
    intensities to set to the \code{trim} and 1-\code{trim} percentile
    intensities when \code{method="winsorize"}.  Default value is 0.05.}
}
\details{
  To summarise the raw data using the default method used by Illumina
  (\code{method="illumina"}) we first remove outliers for each bead
  type on each array. Outliers are beads which have an intensity greater
  than 3 median absolute deviations (MADs) from the bead median intensity on
  the original (un-logged) scale. The \code{n} argument can be changed to
  remove beads with intensity \code{n} MADs above or below the median.
  With outliers removed, the average (mean) intensities of the remaining
  beads are calculated along with the standard error and number of beads.
  
  Other summarisation options are also available. When
  \code{method="mean"}, the average and standard error of all beads
  for a given bead type is calculated on each array.  This would be 
  appropriate if the scanner has been set up to exclude outlier beads  
  from the bead level txt or csv files.
  
  \code{method="median"}, the middle value is returned along with the
  median absolute deviation (rather than standard error) for each bead type.
  When \code{method="trim"}, the trimmed mean and standard error are
  calculated and for \code{method="winsorize"} the winsorised mean and
  standard error are returned.
  
  By setting the \code{log} argument to \code{TRUE}, we calculate outliers and
  summary values on the log2-scale.  

  If there are any NAs or Inf values, they are ignored.
  
  Objects which are created separately by 'createBeadSummaryData' may
  be joined using the \code{combine} function.
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
#produce bead summaries for each array
data(BLData)
BSData = createBeadSummaryData(BLData, log=TRUE, what="G")
dim(BSData)
}
\keyword{methods}

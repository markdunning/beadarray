\name{imageplot}
\alias{imageplot}
\title{imageplot for BeadLevelList object}
\description{
Generates an image plot for data from a \code{BeadLevelList} object.
}
\usage{
imageplot(BLData, array = 1, nrow = 100, ncol = 100, low= NULL,
           high = NULL, ncolors = 123, whatToPlot ="G", log=TRUE,
           zlim=NULL, main=whatToPlot, method="illumina", n = 3,
           trim=0.05, legend=TRUE, ...)
}

\arguments{
  \item{BLData}{\code{BeadLevelList}}
  \item{array}{integer specifying which strip/array to plot}
  \item{nrow}{integer specifying the number of rows to divide the strip/array into}
  \item{ncol}{integer specifying the number of columns to divide the strip/array into}
  \item{low}{colour to use for lowest intensity}
  \item{high}{colour to use for highest intensity}
  \item{ncolors}{The number of colour graduations between high and low}
  \item{whatToPlot}{character string specifying which intensities/values
    to plot. See \code{getArrayData} for a list of possibilities}
  \item{log}{if TRUE, log2 intensities are plotted}
  \item{zlim}{numerical vector of length 2 giving the extreme values of 'z'
    to associate with colours 'low' and 'high'.}
  \item{main}{character string for plot title}
  \item{method}{character string specifying the summarisation method to
    use. Only applicable when \code{whatToPlot="residG"}, \code{"residR"}
    or \code{"residM"}. Refer to the \code{createBeadSummaryData} help
    page for further information.}
  \item{n}{numeric value specifying the number of median absolute deviations (MADs)
    from the median to use as a cut-off for outliers.  The default value is 3. 
    Only applicable when \code{whatToPlot="residG"}, \code{"residR"} or \code{"residM"} 
    and \code{method="illumina"}.  Refer to \code{createBeadSummaryData} help page 
    for further information.}
  \item{trim}{fraction of intensities to remove from the bead summary
    calculations.  Only applicable when \code{whatToPlot="residG"},
      \code{"residR"} or \code{"residM"}.  Refer to
      \code{createBeadSummaryData} help page for further information.}
  \item{legend}{logical, if \code{TRUE}, \code{zlim} and range of data
    is added to plot.}
  \item{\dots}{other graphical parameters to plot that can be specified}
      }
\details{
  Because of the large number of beads on each strip/array, this function
  works by mapping a grid of size specified by the \code{nrow} and \code{ncol}
  arguments and averaging the intensities of the beads within each
  section of the grid.

  The number of rows and columns may change the appearance of the
  plots. If the array is divided into too many squares it will be difficult
  to detect changes. We recommend using \code{nrow=20} and
  \code{ncol=200} for the strips on a BeadChip, and \code{nrow=100}
  \code{ncol=100} for arrays on a SAM.

  An imageplot of the log base 2 foreground intensities is produced by default.
  Other values can be plotted by changing the \code{whatToPlot} argument.
  The default colour scheme ranges from white for low values to blue for high values.
}
\value{

A plot is produced on the current graphical device.
}

\author{Mike Smith, Mark Dunning}
\examples{
data(BLData)
imageplot(BLData)
}
\keyword{hplot}

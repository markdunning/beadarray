\name{imageplot}
\alias{imageplot}
\title{imageplot for BeadArrays object}
\description{
Plots an image plot of a bead level data object.
}
\usage{
imageplot(BLData, array = 1, nrow = 100, ncol = 100, low= NULL,
           high = NULL, ncolors = 123, whatToPlot ="G", log=TRUE,
           zlim=NULL, main=whatToPlot, method="illumina", n = 3,
           trim=0.05, legend=TRUE, ...)
}

\arguments{
  \item{BLData}{Object of class BLData}
  \item{array}{array to plot}
  \item{nrow}{How many rows should the array be divided into}
  \item{ncol}{How many columns should the array be divided into}
  \item{low}{Color for area of lowest intensity}
  \item{high}{color for area of highest intensity}
  \item{ncolors}{The number of color graduations between high and low}
  \item{whatToPlot}{character string specifying which intensities/values
    to plot. See \code{getArrayData} for a list of the possibilities}
  \item{log}{if TRUE, log2 intensities are plotted}
  \item{zlim}{numerical vector of length 2 giving the extreme values of 'z'
    to associate with  colors 'low' and 'high'.}
  \item{main}{character string for plot title}
  \item{method}{chracter string specifying the summarisation method to
    use. Only applicable when \code{whatToPlot="residG"}, \code{"residR"}
    or \code{"residM"}. Refer to help \code{createBeadSummaryData} help
    page for further information.}
  \item{n}{use n median absolute deviations from the median as a
    cut-off for outliers.  Only applicable when
    \code{whatToPlot="residG"}, \code{"residR"} or \code{"residM"}.
    Refer to help \code{createBeadSummaryData} help page for further
      information.}
  \item{trim}{fraction of intensities to remove from the bead summary
    calculations.  Only applicable when \code{whatToPlot="residG"},
      \code{"residR"} or \code{"residM"}.  Refer to
      \code{createBeadSummaryData} help page for further information.}
  \item{legend}{logical, if \code{TRUE}, \code{zlim} and range of data
    is added to plot.}
  \item{\dots}{other graphical parameters to plot that can be specified}
      }
\details{
  Because of the high volume of bead level data, this function
  works by mapping a grid (who's size is specified by the nrow and ncol)
  arguments and averages the intensities of the beads within each
  section of this grid.

  The number of rows and columns may change the appearance of the
  plots. If the array is divided into too many squares it will be difficult
  to detect changes. We recommend using \code{nrow=20} and
  \code{ncol=200} for whole genome arrays, and \code{nrow=100}
  \code{ncol=100} for arrays on a SAM.

  An imageplot of the log base 2 foreround intensities is produced by default.
  Other values can be plotted by changing the \code{whatToPlot} argument.
}
\value{
Plot is produced on current graphical device.
}

\author{Mike Smith, Mark Dunning}
\examples{
data(BLData)
imageplot(BLData)
}
\keyword{hplot}

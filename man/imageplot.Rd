\name{imageplot}
\alias{imageplot}
\title{imageplot for BeadArrays object}
\description{
Plots an image plot of a bead level data object.
}
\usage{
imageplot(BLData, array = 1, nrow = 18, ncol = 2, low= NULL,
           high = NULL, ncolors = 123, whatToPlot ="G", log=TRUE,
           n=3, zlim=NULL, main=whatToPlot, ...)
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
    to plot. Possibilities are \code{"G"}, \code{"Gb"}, \code{"residG"},
    for single channel data and \code{"G"}, \code{"Gb"}, \code{"residG"},
    \code{"R"}, \code{"Rb"}, \code{"residR"}, \code{"M"} (log-ratio) or
    \code{A} (average intensity) for two-colour data}
  \item{log}{if TRUE, log2 intensities are plotted}
  \item{n}{numeric value defining a cut-off of 'n' median absolute
    deviations from the median to determine outliers in
    \code{createBeadSummaryData}. The default value is 3.
    This argument is used to calculate the summarised R or G values
    which are used to calculate the residuals when
    \code{whatToPlot="residR"} or \code{"residG"}.}
  \item{zlim}{numerical vector of length 2 giving the extreme values of 'z'
    to associate with  colors 'low' and 'high'.}
  \item{main}{character string for plot title}
  \item{\dots}{other graphical parameters to plot that can be specified}}
\details{
  Because of the high volume of bead level data, this function
  works by mapping a grid (who's size is specified by the nrow and ncol)
  arguments and averages the intensities of the beads within each
  section of this grid.

  The number of rows and columns may change the appearance of the
  plots. If the array is divided into too many squares it will be difficult
  to detect changes. We recommend using \code{nrow=20} and \code{ncol=200} for whole
  genome arrays, and \code{nrow=100} \code{ncol=100} for arrays on a SAM.

  An imageplot of the log base 2 forgeround intensities is produced by default.
  Other values can be plotted by changing the \code{whatToPlot} argument.
}
\value{
Plot is produced on current graphical device.
}

\author{Mike Smith, Mark Dunning}
\examples{
data(BLData)

imageplot(BLData, nrow=100, ncol=100)
}
\keyword{hplot}

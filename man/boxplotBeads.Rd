\name{boxplotBeads}
\alias{boxplotBeads}
\title{Box plot of bead intensities}
\description{
  Function to produce box plots of the bead intensities from 
  different bead arrays
}
\usage{
boxplotBeads(BLData, whatToPlot = "G", arrays = NULL, log = TRUE,
          varwidth = TRUE, method = "illumina", n = 3, trim = 0.05, ...)
}
\arguments{
  \item{BLData}{a BeadLevelList object}
  \item{whatToPlot}{character string specifying which intensities to plot.
    See \code{getArrayData} for a list of the possibilities.}
  \item{arrays}{numeric vector of arrays to plot}
  \item{log}{if TRUE log2 intensities are plotted}
  \item{varwidth}{logial, indicating whether box widths should be
    proportional to the number of values in each box plot}
  \item{method}{chracter string specifying the summarisation method to
    use (only applicable when \code{whatToPlot="residG"}, \code{"residR"}
    or \code{"residM"}). Refer to help \code{createBeadSummaryData} help
    page for further information.}
  \item{n}{use n median absolute deviations from the median as a
    cut-off for outliers.  Only applicable when
    \code{whatToPlot="residG"}, \code{"residR"} or \code{"residM"}.
    Refer to help \code{createBeadSummaryData} help page for further
      information.}
  \item{trim}{fraction of intensities to remove from the bead summary
    calculations. Only applicable when \code{whatToPlot="residG"},
      \code{"residR"} or \code{"residM"}.  Refer to
      \code{createBeadSummaryData} help page for further information.}
  \item{\dots}{further graphical parameters to the \code{boxplot}
	function from the graphics package}
}
\details{
  Produces box plots of the raw intensities from a bead array
}
\value{
Plot is produced on current graphical device
}

\author{Matt Ritchie}

\examples{
data(BLData)

boxplotBeads(BLData)
}
\keyword{hplot}

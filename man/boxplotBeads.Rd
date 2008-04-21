\name{boxplotBeads}
\alias{boxplotBeads}
\title{Box plot of bead intensities}
\description{
  Function to produce box plots of the bead intensities from 
  selected strips/arrays from a \code{BeadLevelList} object.
}
\usage{
boxplotBeads(BLData, whatToPlot = "G", arrays = NULL, log = TRUE,
          varwidth = TRUE, method = "illumina", n = 3, trim = 0.05, ...)
}
\arguments{
  \item{BLData}{\code{BeadLevelList}}
  \item{whatToPlot}{character string specifying which intensities to plot.
    See \code{getArrayData} for a list of the possibilities.}
  \item{arrays}{integer (scalar or vector) specifying the strips/arrays to plot.
    If \code{NULL}, all the strips/arrays are plotted.}
  \item{log}{if TRUE log2 intensities are plotted}
  \item{varwidth}{logical, indicating whether box widths should be
    proportional to the number of values in each box plot}
  \item{method}{character string specifying the summarisation method to
    use (only applicable when \code{whatToPlot="residG"}, \code{"residR"}
    or \code{"residM"}). Refer to help \code{createBeadSummaryData} help
    page for further information.}
  \item{n}{numeric value defining a cut-off for the number of 
    median absolute deviations (MADs) from the median to use for 
    determining outliers.  The default value is 3.  Only applicable when
    \code{whatToPlot="residG"}, \code{"residR"} or \code{"residM"} and
    \code{method="illumina"}.  Refer to help \code{createBeadSummaryData} 
    help page for further information.}
  \item{trim}{fraction of intensities to remove from the bead summary
    calculations. Only applicable when \code{whatToPlot="residG"},
      \code{"residR"} or \code{"residM"}.  Refer to
      \code{createBeadSummaryData} help page for further information.}
  \item{\dots}{further graphical parameters to the \code{boxplot}
	function from the graphics package}
}
\details{
  Produces box plots of the specified intensities for selected strips/arrays.
}
\value{
A plot is produced on the current graphical device
}

\author{Matt Ritchie}

\examples{
data(BLData)

boxplotBeads(BLData)
}
\keyword{hplot}

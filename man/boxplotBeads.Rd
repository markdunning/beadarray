\name{boxplotBeads}
\alias{boxplotBeads}
\title{Box plot of bead intensities}
\description{
  Function to produce box plots of the bead intensities from 
  different bead arrays
}
\usage{
boxplotBeads(BLData, whatToPlot = "G", arrays = NULL, log = TRUE,
          n=3, varwidth = TRUE, ...)
}
\arguments{
  \item{BLData}{a BeadLevelList object}
  \item{whatToPlot}{character string specifying which intensities to plot.
    Possibilities are \code{"G"}, \code{"Gb"} for single channel data
    and \code{"G"}, \code{"Gb"}, \code{"R"} and \code{"Rb"} for
    two-colour data}
  \item{arrays}{numeric vector of arrays to plot}
  \item{log}{if TRUE log2 intensities are plotted}
  \item{varwidth}{logial, indicating whether box widths should be
    proportional to the number of values in each box plot}
  \item{n}{numeric value defining a cut-off of 'n' median absolute
    deviations from the median to determine outliers in
    \code{createBeadSummaryData}. The default value is 3.
    This argument is used to calculate the summariesed R or G values
    which are used to calculate the residuals when
    \code{whatToPlot="residR"} or \code{"residG"}.}
  \item{\dots}{further graphical parameters to \code{boxplot}
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

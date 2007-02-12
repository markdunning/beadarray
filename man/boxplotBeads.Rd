\name{boxplotBeads}
\alias{boxplotBeads}
\title{Box plot of bead intensities}
\description{
  Function to produce box plots of the bead intensities from 
  different bead arrays
}
\usage{
boxplotBeads(BLData, whatToPlot = "G", arrays = NULL, log = TRUE,
          varwidth = TRUE, ...)
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
  \item{\dots}{further graphical parameters to \code{boxplot}
    function from graphics package}
}
\details{
  Produces box plots or the raw intensities from a bead array}
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

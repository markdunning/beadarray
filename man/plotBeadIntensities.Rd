\name{plotBeadIntensities}
\alias{plotBeadIntensities}
\title{Plot bead intensities}
\description{
  Function to plot the intensities of all beads of a particular
  type on a strip/array.
}
\usage{
plotBeadIntensities(BLData, ProbeIDs, arrays, log = FALSE,
          whatToPlot="G", ProbeCols=NULL, ylim=NULL,...)
}
\arguments{
  \item{BLData}{\code{BeadLevelList}}
  \item{ProbeIDs}{numeric value giving the ProbeID(s) for the bead
    type(s) to plot}
  \item{arrays}{integer (scalar or vector) specifying the strip(s)/array(s) to plot}
  \item{log}{if TRUE log2 intensities are plotted}
  \item{whatToPlot}{character string specifying which intensities to plot.
    Possibilities are \code{"G"}, \code{"Gb"} for single channel data
    and \code{"G"}, \code{"Gb"}, \code{"R"} and \code{"Rb"} for
    two-colour data}
  \item{ProbeCols}{vector of colour names to assign to each bead type}
  \item{ylim}{vector providing the upper and lower bounds for the
    y-axis in the plotting window}
  \item{\dots}{further graphical parameters to \code{plot}}
}
\details{
Boxplots of the specified intensities on the original or log-scale for each bead type specified by \code{ProbeIDs} on selected \code{arrays} is produced.

A vertical line separates the intensities from different strips/arrays.
}

\value{
A plot is displayed on the current graphical device.
}

\author{Mark Dunning}

\examples{
data(BLData)

plotBeadIntensities(BLData, arrays=1:4, ProbeIDs=c(2,1000), log=TRUE)
}
\keyword{hplot}

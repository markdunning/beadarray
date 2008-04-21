\name{plotBeadDensities}
\alias{plotBeadDensities}
\title{plot densities of bead intensities}
\description{
  Function to produce smoothed density plots of the bead intensities from 
  different bead strips/arrays.
}
\usage{
plotBeadDensities(BLData, whatToPlot = "G", arrays = NULL, log = TRUE,
          type="l", col=1, xlab="Intensity", ylab="Density",
          xlim=NULL, ylim=NULL, ...)
}
\arguments{
  \item{BLData}{\code{BeadLevelList}}
  \item{whatToPlot}{character string specifying which intensities 
    to plot.
    Possibilities are \code{"G"}, \code{"Gb"} for single channel data
    and \code{"G"}, \code{"Gb"}, \code{"R"} and \code{"Rb"} for
    two-colour data}
  \item{arrays}{integer (scalar or vector) specifying the strip(s)/array(s) to plot. 
    If set to \code{NULL} (default value), data from all strips/arrays 
    are plotted}
  \item{log}{if TRUE log2 intensities are plotted}
  \item{type}{character string specifying what type of plot to draw}
  \item{col}{the colours for lines and points}
  \item{xlab}{label for x-axis}
  \item{ylab}{label for y-axis}
  \item{xlim}{numeric vector specifying x-axis limits}
  \item{ylim}{numeric vector specifying y-axis limits}
  \item{\dots}{further graphical parameters to \code{plot}}
}
\details{
  Produces density plots of the raw intensities from a \code{BeadLevelList}
  When \code{arrays="NULL"}, densities from all arrays are plotted.
}
\value{
A plot is produced on the current graphical device
}

\author{Matt Ritchie}

\examples{
data(BLData)
plotBeadDensities(BLData)
}
\keyword{hplot}

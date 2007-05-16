\name{plotBeadDensities}
\alias{plotBeadDensities}
\title{plot density of bead intensities}
\description{
  Function to produce smoothed density plots of the bead intensities from 
  different bead arrays
}
\usage{
plotBeadDensities(BLData, whatToPlot = "G", arrays = NULL, log = TRUE,
          type="l", col=1, xlab="Intensity", ylab="Density",
          xlim=NULL, ylim=NULL, ...)
}
\arguments{
  \item{BLData}{a BeadLevelList object}
  \item{whatToPlot}{character string specifying which intensities 
    to plot.
    Possibilities are \code{"G"}, \code{"Gb"} for single channel data
    and \code{"G"}, \code{"Gb"}, \code{"R"} and \code{"Rb"} for
    two-colour data}
  \item{arrays}{numeric vector of arrays to plot. If set to 
    \code{NULL} (default value), data from all arrays are plotted}
  \item{log}{if TRUE log2 intensities are plotted}
  \item{type}{character string specifying what type of plot to draw}
  \item{col}{the colors for lines and points}
  \item{xlab}{label for x-axis}
  \item{ylab}{label for y-axis}
  \item{xlim}{numeric vector specifying x-axis limits}
  \item{ylim}{numeric vector specifying y-axis limits}
  \item{\dots}{further graphical parameters to \code{plot}}
}
\details{
  Produces density plots of the raw intensities from a BeadArray
  When \code{arrays="NULL"}, densities from all arrays are plotted.
}
\value{
Plot is produced on current graphical device
}

\author{Matt Ritchie}

\examples{
data(BLData)
plotBeadDensities(BLData)
}
\keyword{hplot}

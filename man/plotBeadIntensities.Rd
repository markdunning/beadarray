\name{plotBeadIntensities}
\alias{plotBeadIntensities}
\title{Plot bead intensities}
\description{
  Function to plot the intensities of all beads of a particular
  type on an array.
}
\usage{
plotBeadIntensities(BLData, ProbeIDs, arrays, log = FALSE,
          whatToPlot="G", ProbeCols=NULL, ylim=NULL,...)
}
\arguments{
  \item{BLData}{a BeadLevelList object}
  \item{ProbeIDs}{numeric value giving the ProbeID(s) for the bead
    type(s) we are interested in}
  \item{arrays}{numeric value giving the array we are interested in}
  \item{log}{if TRUE log2 intensities are plotted (and used to
    find outliers)}
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
Boxplots of the G, Gb, R or Rb intensities on the original or log-scale for each bead type specified by \code{ProbeIDs} on selected \code{arrays} is produced.
}
\value{
Plot is produced on current graphical device
}

\author{Mark Dunning}

\examples{
data(BLData)

plotBeadIntensities(BLData, arrays=seq(1:4), ProbeIDs=c(2,1000), log=TRUE)
}
\keyword{hplot}

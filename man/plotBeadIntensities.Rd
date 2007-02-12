\name{plotBeadIntensities}
\alias{plotBeadIntensities}
\title{Plot bead intensities}
\description{
  Function to plot the intensities of all beads of a particular
  type on an array against their distance from the centre of the
  array. Any outliers for the bead type will be highlighted.
}
\usage{
plotBeadIntensities(BLData, ProbeIDs, arrays, log = FALSE, n=3,
          whatToPlot="G", ProbeCols=NULL, ylim=NULL,...)
}
\arguments{
  \item{BLData}{a BeadLevelList object}
  \item{ProbeIDs}{numeric value giving the ProbeID for the bead
    type we are interested in}
  \item{arrays}{numeric value giving the array we are interested in}
  \item{log}{if TRUE log2 intensities are plotted (and used to
    find outliers)}
  \item{n}{outlier values for the bead type are beads with intensity
    more than 'n' times away from the mean}
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
We find the intensities for each bead of the specified bead type by
finding the G, Gb, R or Rb intensities for each bead in the array
whose BLData\$ProbeID matches the specified ProbeID. The distance from
the centre of the hexagon can be found in a similar fashion by looking
at the GrnX and GrnY values.

The outliers for the bead type are beads which have intensity more
than 'n' times away from the standard deviation using logged or
un-logged bead intensities (specified using 'log' argument).

The mean for the bead type and lines 3 standard deviations from
the mean are shown on the plot.
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
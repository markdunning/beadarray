\name{plotProbeDistances}
\alias{plotProbeDistances}

\title{Plot the distances between a set of beads}

\description{
The function calculates the distance matrix for a set of bead co-ordinates and plots in the form of a histogram. 
}

\usage{
plotProbeDistances(BLData, probe, array)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{BLData}{BeadLevelList object containing bead-level data}
\item{probe}{numeric value specifying the ProbeID}
\item{array}{numeric value specifying the array}

}

\details{

A histogram is produced with the mean of the distance matrix indicated in red. A normal distribution with mean = 650 and sd =350 is overlayed.

}

\value{
Plot is produced on the current graphical device
}


\author{Mark Dunning}

\examples{

data(BLData)
plotProbeDistances(BLData, 2, 1)

}
\seealso{
\code{\link{getProbeCoords}}
\code{\link{plotDistances}}
}
\keyword{methods}
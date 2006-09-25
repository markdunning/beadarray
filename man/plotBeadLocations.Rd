\name{plotBeadLocations}
\alias{plotBeadLocations}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Show bead locations}
\description{
Shows location of a set of beads on an array. The beads can either be defined to be all beads with particular ProbeIDs or individual beadIDs
}
\usage{
plotBeadLocations(BLData, ProbeIDs=NULL, BeadIDs=NULL, array, SAM=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{a BeadLevelList object}
  \item{ProbeIDs}{a vector of ProbeIDs to plot}
  \item{BeadIDs}{a vector of beadIDs to plot}
  \item{array}{numeric value specifying the array}
  \item{SAM}{if TRUE then the data is assumed to be taken from a SAM array and therefore hexagonal}
  \item{\dots}{further graphical parameters to \code{plot}}
}
\details{

The outline of the hexagonal array is drawn and the locations of the specified beads are overlayed. 

}
\value{
Plot is produced on the current graphical device.
}

\examples{

data(BLData)

#plot all beads with ProbeID 2 on array 1
plotBeadLocations(BLData, ProbeIDs=2, array=1)

#plot all beads with either ProbeID 2 or 10

plotBeadLocations(BLData, ProbeIDs=c(2, 10), array=1)

#can also change plotting characters

plotBeadLocations(BLData, ProbeIDs=c(2, 10), array=1, col="red", pch="X", cex=2)

#find all outliers on the first array and plot their locations
o=findAllOutliers(BLData, 1)

plotBeadLocations(BLData, beadIDs=o, array=1)



}

\author{Mark Dunning}


\keyword{hplot}

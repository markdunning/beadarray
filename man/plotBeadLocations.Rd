\name{plotBeadLocations}
\alias{plotBeadLocations}
\alias{plotCoords}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Show bead locations}
\description{
Shows location of a set of beads on an array. The beads can either be defined to be all beads with particular probeIDs or individual beadIDs
}
\usage{
plotBeadLocations(BLData, probeIDs=NULL, beadIDs=NULL, array, label=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{a BeadLevelList object}
  \item{probeIDs}{a vector of probeIDs to plot}
  \item{beadIDs}{a vector of beadIDs to plot}
  \item{array}{numeric value specifying the array}
  \item{label}{if TRUE every bead plotted with by given a numeric identifier on the plot}		
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

#plot all beads with probeID 2 on array 1
plotBeadLocations(BLData, probeIDs=2, array=1)

#plot all beads with either probeID 2 or 10

plotBeadLocations(BLData, probeIDs=c(2, 10), array=1)

#can also change plotting characters

plotBeadLocations(BLData, probeIDs=c(2, 10), array=1, col="red", pch="X", cex=2)

#find all outliers on the first array and plot their locations
o=findAllOutliers(BLData, 1)

plotBeadLocations(BLData, beadIDs=o, array=1)



}

\author{Mark Dunning}

\seealso{\code{\link{findAllOutliers}}}
\keyword{hplot}

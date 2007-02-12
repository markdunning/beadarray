\name{plotBeadLocations}
\alias{plotBeadLocations}
\title{Show bead locations}
\description{
  Shows location of a set of beads on an array. The beads can
  either be defined to be all beads with particular ProbeIDs
  or individual beadIDs
}
\usage{
plotBeadLocations(BLData, ProbeIDs = NULL, BeadIDs = NULL, array = 1,
                   SAM = FALSE, xlab = "x-coordinate",
                   ylab = "y-coordinate",
                   main = paste("Bead", ProbeIDs,"locations"), ...)
}

\arguments{
  \item{BLData}{a BeadLevelList object}
  \item{ProbeIDs}{a vector of ProbeIDs to plot}
  \item{BeadIDs}{a vector of beadIDs to plot}
  \item{array}{numeric value specifying the array}
  \item{SAM}{if TRUE then the data is assumed to be
    taken from a SAM array and therefore hexagonal}
  \item{xlab}{character string specifying x-axis label}
  \item{ylab}{character string specifying y-axis label}
  \item{main}{character string specifying plot title}
  \item{\dots}{further graphical parameters to \code{plot}}
}
\details{
  The outline of the hexagonal array is drawn and the
  locations of the specified beads are overlayed. 
}
\value{
Plot is produced on the current graphical device.
}

\examples{
data(BLData)

#plot all beads with ProbeID 2 on array 1
plotBeadLocations(BLData, array=1, ProbeIDs=2, SAM=TRUE)

#find all outliers on the first array and plot their locations
o=findAllOutliers(BLData, 1)
plotBeadLocations(BLData, BeadIDs=o, array=1, SAM=TRUE)
}

\author{Mark Dunning}

\keyword{hplot}
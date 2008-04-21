\name{plotBeadLocations}
\alias{plotBeadLocations}
\title{Plot bead locations}
\description{
  Shows location of a set of beads on a strip/array. The beads can
  either be defined to be all beads with particular ProbeIDs
  or as rows in \code{BeadLevelList}.
}
\usage{
plotBeadLocations(BLData, ProbeIDs = NULL, BeadIDs = NULL, array = 1,
                   SAM = FALSE, xlab = "x-coordinate",
                   ylab = "y-coordinate",
                   main = paste("Bead", ProbeIDs,"locations"), ...)
}

\arguments{
  \item{BLData}{\code{BeadLevelList}}
  \item{ProbeIDs}{a vector of ProbeIDs to plot}
  \item{BeadIDs}{logical/integer vector specifying which rows of data from
   \code{BeadLevelList}) to plot (used if \code{ProbeIDs} is \code{NULL})}
  \item{array}{integer specifying the strip/array to plot}
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
A plot is produced on the current graphical device.
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
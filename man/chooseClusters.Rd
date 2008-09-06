\name{chooseClusters}
\alias{chooseClusters}
\title{Choose Clusters}
\description{
Find large clusters of beads. 
}
\usage{
chooseClusters(IDs, neighbours, cutoff = 8)
}

\arguments{
  \item{IDs}{IDs of beads to be clustered.}
  \item{neighbours}{A Neighbours matrix - obtained from \code{\link{generateNeighbours}}.}
  \item{cutoff}{Integer - threshold for the minimum size a cluster must be.}
}

\details{
	This function will find which beads are in large clusters. Using a flood fill algorithm, it finds clusters of beads, determines the size of each, and then returns only the beads in clusters of size greater than \code{cutoff}. It is primarily used in \code{\link{BASHCompact}} and \code{\link{BASHDiffuse}}.
}

\value{
Vector of bead IDs. (This will be a subset of the argument \code{IDs})
}

\author{Jonathan Cairns}

\seealso{\code{\link{BASHCompact}}, \code{\link{BASHDiffuse}}, \code{\link{closeImage}}}

\examples{
data(BLData)
neighbours <- generateNeighbours(BLData,1)
o <- findAllOutliers(BLData,1,log = TRUE)
clusters8 <- chooseClusters(o, neighbours)
clusters12 <- chooseClusters(o, neighbours, cutoff = 12) ## only larger clusters

x11()
plotBeadLocations(BLData,array=1,BeadIDs = clusters8, pch = ".")
}
\keyword{misc}

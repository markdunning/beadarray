\name{closeImage}
\alias{closeImage}
\title{Close Image}
\description{
Find the closure of a set of beads on an array.
}
\usage{
closeImage(IDs, neighbours, cinvasions = 10)
}

\arguments{
  \item{IDs}{IDs of beads to be closed.}
  \item{neighbours}{A Neighbours matrix - obtained from \code{\link{generateNeighbours}}.}
  \item{cinvasions}{The number of invasions used when dilating and eroding.}
      }

\details{
	This function "closes" the set of beads supplied, as used in the \code{\link{BASH}} functions. It dilates (expands) the image, and then erodes (contracts) it. Each is done via an invasion process - if we let the set of beads supplied be called S, then dilation considers all neighbours of beads in S, and adds them to S. Erosion finds all beads in S with neighbours outside of S, and removes them from S.

The result of this process is to close "holes" in the group of specified beads during the dilation. These are not reopened during the erosion.
}

\value{
An updated vector of bead IDs (of which the argument \code{IDs} will be a subset).
}

\author{Jonathan Cairns}

\seealso{\code{\link{generateNeighbours}}}

\examples{
data(BLData)

##This process is equivalent to one iteration of BASHCompact.
##o <- findAllOutliers(BLData,4)
##neighbours <- generateNeighbours(BLData,4)
##o.clusters <- chooseClusters(o, neighbours)
##o.compact <- closeImage(o.clusters, neighbours)
}

\keyword{misc}

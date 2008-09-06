\name{generateNeighbours}
\alias{generateNeighbours}
\title{Generate Neighbours Matrix for BeadLevelList object}
\description{
Generates a Neighbours matrix from the X and Y co-ordinates in a \code{BeadLevelList} object.
}
\usage{
generateNeighbours(BLData, array, thresh = 2.2, margin = 10, window = 30)
}

\arguments{
  \item{BLData}{\code{BeadLevelList}}
  \item{array}{integer specifying which strip/array to plot}
  \item{window}{numeric value, specifying window size (see below)}
  \item{margin}{numeric value, specifying size of window margin (see below)}
  \item{thresh}{numeric value, which determines how large links are removed. (see below)}
      }

\details{
	generateNeighbours determines, for each bead on the array, which beads are next to it. It assumes that the beads are in a hexagonal lattice.

	The algorithm used first links each bead to its 6 closest neighbours. It then removes the longest link if its squared length is more than \code{thresh} multiplied by the squared length of the next longest link. A similar process is applied to the 2nd and 3rd longest links.

	Finally, any one way links are removed (i.e. a link between two beads is only preserved if each bead considers the other to be its neighbour).

	To ease computation, the algorithm only computes neighbours of beads in a square window of side length \code{2*(window)} which travels across the array. Beads in a margin around the square, of width \code{(margin)}, are also considered as possible neighbours.

	The Neighbours matrix is designed for use with the BASH functions.
}

\value{
A matrix with 6 columns, and a number of rows equal to the number of beads on the array. The neighbours of bead i are found in row i. 0 denotes a deleted link. (For example, if row 15 consists of 20, 35, 0, 0, 60, 4, then bead 15 is considered a neighbour of beads 4, 20, 35 and 60.)
}

\author{Jonathan Cairns}
\examples{
data(BLData)
neighbours <- generateNeighbours(BLData, 1)
}
\keyword{misc}

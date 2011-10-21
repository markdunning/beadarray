\name{generateNeighbours}
\alias{generateNeighbours}
\title{Generate matrix of neighbouring beads}
\description{
Generates a neighbours matrix from either a .locs file or the X and Y coordinates in a \code{beadLevelData} object.
}
\usage{
generateNeighbours(BLData, array = 1, useLocs = TRUE, window = 30, margin = 10, thresh = 2.2)
}

\arguments{
  \item{BLData}{An object of class \code{\link{beadLevelData-class}}}
  \item{array}{integer specifying which section/array to process}
  \item{useLocs}{logical value, specifiying whether the .locs file (if present) should be used to determine neighbours.}
  \item{window}{numeric value, specifying window size (see below)}
  \item{margin}{numeric value, specifying size of window margin (see below)}
  \item{thresh}{numeric value, which determines how large links are removed. (see below)}
}

\details{
generateNeighbours determines, for each bead on the array, which beads are next to it. It assumes that the beads are in a hexagonal lattice.

If the \code{.locs} file is present and \code{useLocs = TRUE} then the ordering of the \code{.locs} file is used to infere the grid layout.  This is far quicker than the alternative and is thus recommended, but can only be used on BeadChip platforms.  If the data is from a Sentrix Array \code{useLocs} is automatically set to \code{FALSE} and the following algorithm is applied instead.

The algorithm used first links each bead to its 6 closest neighbours. It then removes the longest link if its squared length is more than \code{thresh} multiplied by the squared length of the next longest link. A similar process is applied to the 2nd and 3rd longest links.

Finally, any one way links are removed (i.e. a link between two beads is only preserved if each bead considers the other to be its neighbour).

To ease computation, the algorithm only computes neighbours of beads in a square window of side length \code{2*(window)} which travels across the array. Beads in a margin around the square, of width \code{(margin)}, are also considered as possible neighbours.
}

\value{
A matrix with 6 columns.  Each row corresponds to a bead in the passed \code{\link{beadLevelData-class}} and the six entries are the indices of the 6 neighbouring beads.  Values of NA indicate that the neighbouring bead appears to be missing, either due to failing Illumina's decoding or being at the edge of the array.
}

\author{Jonathan Cairns, Mike Smith}

\seealso{\code{\link{HULK}}, \code{\link{BASH}}}

\examples{

\dontrun{
if(require(beadarrayExampleData)){

data(exampleBLData);
neighbours <- generateNeighbours(exampleBLData, array = 1, useLocs = FALSE);
}

}

}

\references{
Lynch AG, Smith ML, Dunning MJ, Cairns JM, Barbosa-Morais NL, Tavare S. 
beadarray, BASH and HULK - tools to increase the value of Illumina BeadArray experiments. 
In A. Gusnato, K.V. Mardia, & C.J. Fallaize (eds), Statistical Tools for Challenges in Bioinformatics. 2009 pp. 33-37. Leeds, Leeds University Press.

}
\keyword{misc}

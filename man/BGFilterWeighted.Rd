\name{BGFilterWeighted}
\alias{BGFilterWeighted}
\title{Weighted Background Filter}
\description{
Finds local weighted means at each bead.
}
\usage{
BGFilter(E = NULL, neighbours, invasions = 20, weights = NULL)
}

\arguments{
  \item{E}{Error Image}
  \item{neighbours}{A Neighbours matrix. Required.}
  \item{invasions}{Integer - Number of invasions. This argument is passed to the function BGfilter.}
  \item{weights}{Numerical vector - A vector of weights, from 0 to 1, to consider in the analysis. (see below.)}
}

\details{
	This function finds the weighted mean of local bead intensities, using intensities from the given error image.

	To obtain our weighted mean for each bead, we use an invasion process. Links between beads are defined in the neighbours matrix. We define the local beads as those which can be reached in \code{invasions} steps from the first bead, and then collect their error values.

	We take a weighted mean of these error values, where the weights are calculated by taking the product of:
	a) 1/(the number of steps required to get to the bead from the central bead)
	b) (if supplied) the weights defined through the \code{weights} parameter.

	This weighted mean is then assigned to the central bead.
}

\value{
A vector - the weighted means. (NB: Whilst \code{BGFilter} manipulates the error image and returns an updated error image, e.g. subtracting the local median, this function does not - it merely returns the local weighted means.)
}

\author{Jonathan Cairns}

\seealso{\code{\link{BGFilter}}}

\examples{
data(BLData)
E <- generateE(BLData,1,method = "mean")
neighbours <- generateNeighbours(BLData,1)
##bgf <- BGFilterWeighted(E, neighbours)
}
\keyword{misc}

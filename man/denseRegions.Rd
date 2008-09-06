\name{denseRegions}
\alias{denseRegions}
\title{Find Dense Regions of Points (as used in diffuse defect analysis).}
\description{
Given a list of beads, this function finds dense regions of beads on the list.
}
\usage{
denseRegions(IDs, neighbours, ignore = NULL, sig = 0.0001, invasions = 10)
}

\arguments{
  \item{IDs}{Vector - IDs of beads to find dense regions of.}
  \item{neighbours}{A Neighbours matrix - obtained from \code{\link{generateNeighbours}}.}
  \item{ignore}{Vector - IDs of beads to be ignored during this process.}
  \item{sig}{Significance of the Binomial test performed within each kernel.}
  \item{invasions}{Integer - No of invasions used to generate the kernel.}
}

\details{
	This function, given a list of bead IDs, finds regions where these marked beads are denser.

	To do this, we use a "sliding kernel" technique. For each bead, we find the "kernel", a local neighbourhood of beads, obtained via invasion along links defined in the neighbours matrix. We count the number of beads in the kernel, and we also count how many of these are beads are marked. Now we test the density of this region with a binomial test.

	Assuming that we expect the marked beads to be evenly distributed across the array, then the number of marked beads in the kernel should have distibution Bin(n,p) under the null hypothesis, where n is the total number of beads in the kernel, and p is the proportion of marked beads on the entire array. We test this hypothesis at a level defined by \code{sig}, and on rejection of the null hypothesis we label the kernel's central bead as being part of a dense region. This is performed for the kernel about each bead.

	If \code{ignore} is specified, then these beads will be completely removed before analysis. Any links attached to a removed bead are severed.
}

\value{
Vector - IDs of beads in dense regions.
}

\author{Jonathan Cairns}

\seealso{\code{\link{generateNeighbours}}, \code{\link{BASHDiffuse}}}


\examples{
data(BLData)
E <- generateE(BLData,1)
E <- generateE(BLData,1, invasions = 10) #reduced no of invasions to increase speed.
E <- generateE(BLData,1, bgfilter = "none") #residuals (median)
}
\keyword{misc}

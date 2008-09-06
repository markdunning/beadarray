\name{BGFilter}
\alias{BGFilter}
\title{Background Filter}
\description{
Performs various image transforms, based on statistics from local beads.
}
\usage{
BGFilter(E = NULL, neighbours, invasions = 20, method = "median")
}

\arguments{
  \item{E}{Error Image}
  \item{neighbours}{A Neighbours matrix. Required.}
  \item{invasions}{Integer - Number of invasions. This argument is passed to the function BGfilter.}
  \item{method}{Method for computing local statistics. Options are \code{"median"}, \code{"mean"}, \code{"MAD"}, \code{"medianMAD"}}
}

\details{
	This function transforms an error image based on a local statistic.

	To obtain our statistic, we use an invasion process. Links between beads are defined in the neighbours matrix. We define the local beads as those which can be reached in \code{invasions} steps from the first bead, and then collect their values.

	\code{method = "median"} subtracts the local median from each error intensity.

	\code{method = "mean"} subtracts the local mean from each error intensity.

	\code{method = "MAD"} divides each bead's intensity by the MAD (median absolute devation from the median) of local beads.

	\code{method = "medianMAD"} subtracts the local median from each error intensity, and then divides each intensity by the local MAD.
}

\value{
A vector - the updated error image.
}

\author{Jonathan Cairns}

\seealso{\code{\link{BGFilter}}}

\examples{
data(BLData)
E <- generateE(BLData,1,bgfilter = "none")
neighbours <- generateNeighbours(BLData,1)
E.MAD <- BGFilter(E, neighbours, method = "MAD")
E.median <- BGFilter(E, neighbours, method = "median")
}
\keyword{misc}

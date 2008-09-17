\name{generateE}
\alias{generateE}
\title{Generate Error Image for BeadLevelList object}
\description{
Generates an Error Image from the data in a \code{BeadLevelList} object.
}
\usage{
generateE(BLData, array, neighbours = NULL, log = TRUE, method = "median", what = "residG", bgfilter = "none", invasions = 20)
}

\arguments{
  \item{BLData}{\code{BeadLevelList}}
  \item{array}{integer specifying which strip/array to plot}
  \item{neighbours}{A Neighbours matrix. Optional - if left NULL, it will be computed.}
  \item{log}{Logical. If TRUE, compute residuals on the log scale.}
  \item{method}{Method for computing residuals. Options are \code{"mean"} and \code{"median"}}
  \item{what}{What to derive the error image from, as used in \code{getArrayData}.}
  \item{bgfilter}{Method passed to the function \code{\link{BGFilter}}. Options are "none", "mean", "median", "MAD" and "medianMAD".}
  \item{invasions}{Integer - Number of invasions. This argument is passed to the function BGfilter.}
}

\details{
	generateE creates an error image, usually based on bead residuals. This output can then be fed into \code{BASHDiffuse} or \code{\link{BASHExtended}}.

	If \code{what} is \code{residG}, \code{residR}, or \code{residM}, then residuals are calculated based on \code{method}. For other values of \code{what}, the residuals are not calculated.

	We then apply a "background filter" to this data, using the function \code{\link{BGFilter}} with arguments \code{bgfilter} and \code{invasions} - see its help file for more details. The background filter subtracts an estimate of the local background of the error image, and/or scales by the local MAD. This step is disabled by using \code{bgfilter = "none"}.
}

\value{
An "Error Image" - a vector of length equal to the number of beads on the array.
}

\author{Jonathan Cairns}

\seealso{\code{\link{BGFilter}}}


\examples{
data(BLData)
E <- generateE(BLData,1)
E <- generateE(BLData,1, invasions = 10) #reduced no of invasions to increase speed.
E <- generateE(BLData,1, bgfilter = "none") #residuals (median)
}
\keyword{misc}

\name{BASHExtended}
\alias{BASHExtended}
\title{BASH - Extended Defect Analysis}
\description{
Returns a score, which assesses the extent to which the background is changing across the array/strip.
}
\usage{
BASHExtended(BLData, array, transFun = logGreenChannelTransform, neighbours = NULL, useLocs = TRUE, E = NULL, E.BG = NULL)
}

\arguments{
  \item{BLData}{\code{BeadLevelList}}
  \item{array}{integer specifying which strip/array to plot}
  \item{transFun}{Function to use to transform data prior to running BASH.}
  \item{neighbours}{A Neighbours matrix. Optional - if left NULL, it will be computed, using default \code{\link{generateNeighbours}} settings.}
  \item{useLocs}{Logical value, specifiying whether the .locs file (if present) should be used to determine neighbours.}
  \item{E}{Numerical vector - The error image to use. Optional - if left blank, it will be computed, using \code{generateE} (with \code{bgfilter = "none"}, i.e. no background filter applied).}
  \item{E.BG}{Numerical vector - The background error image to use. Optional - if left blank, it will be computed from \code{E}, using default \code{BGFilter} settings (i.e. \code{method = "median"}).}
}

\details{
	\code{BASHExtended} assesses the change of background across an array.

	The error image used should not be background filtered (as opposed to the error image used in \code{\link{BASHDiffuse}}). Here, E is the error image
}

\value{
Scalar (Extended defect score)
}

\author{Jonathan Cairns}

\seealso{\code{\link{BASH}}, \code{\link{generateNeighbours}}, }

\examples{
data(BLData)
sn <- sectionNames(BLData)
extended <- NULL

for(i in 1:length(sn))
{
	extended[i] <- BASHExtended(BLData, i)
}
}
\references{
	Mayte Suarez-Farinas, Maurizio Pellegrino, Knut M. Wittkwosky and Marcelo O. Magnasco (2007). Harshlight: A "corrective make-up" program for microarray chips. R package version 1.8.0. http://asterion.rockefeller.edu/Harshlight/

}
\keyword{misc}

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
  \item{array}{integer specifying which section/array to analyse}
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


\dontrun{


if(require(beadarrayExampleData)){


 data(exampleBLData)
 extended <- BASHExtended(exampleBLData, 1)


}

}
}
\references{
	J. M. Cairns, M. J. Dunning, M. E. Ritchie, R. Russell, and A. G. Lynch (2008). BASH: a tool for managing BeadArray spatial artefacts. Bioinformatics 15; 24(24)

}
\keyword{misc}

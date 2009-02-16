\name{HULKresids}
\alias{HULKResids}
\title{HULK - Residuals}
\description{
Calculates an set of weighted average residuals, one for each probe, based on the residuals of the surrounding probes.
}
\usage{
HULKResids(BLData, array, neighbours = NULL, invasions = 20, what = "G")
}

\arguments{
  \item{BLData}{\code{\link{BeadLevelList}}}
  \item{array}{integer specifying which strip/array to plot}
  \item{neighbours}{A Neighbours matrix. Optional - if left NULL, it will be computed.}
  \item{invasions}{Integer - Number of invasions used when identifying
    neighbouring beads.}
  \item{what}{Specify the data in the BLData to create the residuals
  from.  Defaults to the foreground intensities of the green channel.}
}

\details{
  HULKResids calculates a weighted average residual for each probe on
  the specified array of BLData.  It makes use of the same neighbourhood
  calculations as other BASH functions.  The average
  residuals are calculated as the sum of the residuals for each bead in
  the neighbourhood, divided by 1 plus the number of invasions it took to reach
  that bead.  It is intended that \code{HULKResids} be called through
  \code{\link{HULK}}, but it is quite possible to call it as a stand alone function.
}

\value{
  A vector containing an average residual for each bead on the specified
  array of BLData.
}

\author{Mike Smith}

\seealso{\code{\link{HULK}}, \code{\link{BASH}}}

\examples{
data(BLData)
o <- HULKResids(BLData, 1)
}
\references{
Cairns JM, Dunning MJ, Ritchie ME, Russell R, Lynch AG. (2008). BASH: a
tool for managing BeadArray spatial artefacts. Bioinformatics, 24(24):2921-2.

}
\keyword{misc}

\name{HULK}
\alias{HULK}
\title{HULK - Bead Array Normalization by NEighbourhood Residuals}
\description{
Normalizes an probe intensities by calculating a weighted average
residual based on the residuals of the surrounding probes.
}
\usage{
HULK(BLData, array, neighbours = NULL, invasions = 20, what = "G")
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
  HULK is a method of intensity normalization based upon the
  BASH framework.  Firstly For each bead a local neighbourhood of beads is
  determined, using the same process as the other BASH functions.

  For each bead a weighted average residual is calculated.  The average
  residual is calculated as the sum of the residuals for each bead in
  the neighbourhood, divided by 1 plus the number of invasions it took to reach
  that bead.  This calculation is made by a call to
  \code{\link{HULKResids}}.

  The average residuals are then subtracted from each bead and the
  resulting \code{\link{BeadLevelList}} object is returned.
}

\value{
  An object of class \code{BeadLevelList}
}

\author{Mike Smith}

\seealso{\code{\link{HULKResids}}, \code{\link{BASH}}}

\examples{
data(BLData)
o <- HULK(BLData, 1)
}
\references{
Cairns JM, Dunning MJ, Ritchie ME, Russell R, Lynch AG. (2008). BASH: a
tool for managing BeadArray spatial artefacts. Bioinformatics, 24(24):2921-2.

}
\keyword{misc}

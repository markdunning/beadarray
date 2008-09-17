\name{viewBeads}
\alias{viewBeads}
\title{View Beads}
\description{
View an image of the beads in a certain region, optionally with links between neighbours or certain beads highlighted.
}
\usage{
viewBeads(BLData, array, x, y, xwidth = 100, ywidth = 100, neighbours = NULL, mark = NULL, markcol = "blue", markpch = 21, inten = TRUE, low = "black", high = "green", what = "G", log = TRUE, zlim = NULL, ...)
}

\arguments{
  \item{BLData}{\code{BeadLevelList}}
  \item{array}{integer specifying which strip/array to plot}
  \item{x}{numeric value - x co-ordinate to centre on}
  \item{y}{numeric value - y co-ordinate to centre on}
  \item{xwidth}{numeric value - width of square}
  \item{ywidth}{numeric value - width of square}
  \item{neighbours}{Neighbours matrix (optional) - if specified, links will be drawn between neighbours. (See \code{\link{generateNeighbours}}.)}
  \item{mark}{integral vector (optional) - a list of beadIDs to highlight.}
  \item{markcol}{The colour used for the highlighted beads.}
  \item{markpch}{The colour used for the highlighted beads.}
  \item{inten}{logical - if true, plot coloured circles, with shades corresponding to intensities. Intensities are retrieved using \code{\link{getArrayData}}, and the arguments below.}
  \item{low}{Colour used for low intensities.}
  \item{high}{Colour used for high intensities.}
  \item{what}{Data to be used - passed to \code{\link{getArrayData}}.}
  \item{log}{Data to be used - passed to \code{\link{getArrayData}}.}
  \item{zlim}{Limits to be used. Supply in form c(0,5).}
  \item{...}{Additional arguments passed to \code{\link{getArrayData}}.}
  }

\details{
	viewBeads plots the beads within the defined square region.

	Specifying a \code{neighbours} matrix will result in links between neighbours being plotted. Specifying a \code{mark} vector of beadIDs will result in the beads with these beadIDs being highlighted with a blue circular border.
}

\value{
Outputs to the active graphical device.
}

\author{Jonathan Cairns}

\seealso{\code{\link{generateNeighbours}}}

\examples{
##data(BLData)
##o <- findAllOutliers(BLData,2)

##x11()
##viewBeads(BLData, 2, 1000,1200,250,250, mark = o) ## outliers in this ##region are marked in blue
##x11()
##viewBeads(BLData, 2, 1000,1200,250,250, mark = o, inten = FALSE) ## ##removing the intensities makes the outlier locations easier to see

##neighbours <- generateNeighbours(BLData, 2)
##x11()
##viewBeads(BLData, 2, 1000,1200,250,250, neighbours) ## observe that ##there are many missing beads in this region, affecting the neighbours ##algorithm
}

\keyword{hplot}

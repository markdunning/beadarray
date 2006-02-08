\name{getProbeCoords}
\alias{getProbeCoords}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find Bead Co-ordinates}
\description{
  Finds the x and y co-ordinates of every bead of a given bead type on a particular array
}
\usage{
getProbeCoords(BLData, probe, array)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{BLDataLlist object}
  \item{probe}{numeric value for the ProbeID we want to find the co-ordinates for}
  \item{array}{numeric value for the array we want}
}
\details{
Function simply uses the ProbeID to look up the BLData\$x and BLData\$y values for beads whose BLData\$ProbeID match the given 'ProbeID' argument.
}
\value{
  matrix of numeric values with two columns ($x$ and $y$) and number of rows being the number of beads with this ProbeID present in the array. The $x$ and $y$
values returned are standardised by subtracting the minimum $x$ and $y$ values of any bead on the array.
}

\author{Mark Dunning}



\examples{
data(BLData)
getProbeCoords(BLData, 2, 1)

}
\keyword{methods}
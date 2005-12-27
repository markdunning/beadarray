\name{getDimensions}
\alias{getDimensions}

\title{Get dimension of an array}
\description{
Returns the maximum x and y co-ordinates for an array
}

\usage{
getDimensions(BLData, array)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{BLData}{\code{BeadLevelList} object containing bead-level data}
\item{array}{numeric value for the array required}}

\details{

}
\value{
List of length 2 giving the maximum x and y co-ordinates respectively
}


\author{Mark Dunning}

\examples{

data(BLData)
getDimensions(BLData, 1)

}

\keyword{methods}
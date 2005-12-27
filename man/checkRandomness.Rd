\name{checkRandomness}
\alias{checkRandomness}
\alias{chitest}
\title{Test for random distribution of co-ordinates}

\description{
Function tests to see if a set of co-ordinates are randomly positioned around an array. This can be used to test the distribution of beads of the same 
type on the same array or outliers for an array.
}

\usage{
checkRandomness(BLData, array, coords)
}

\arguments{
\item{BLData}{BeadLevelList containing bead-level data}
\item{array}{numeric value specifying the array we want to check the co-ordinates on}
\item{coords}{n by 2 vector where first column is list of x co-ordinates and second column is list of y co-ordinates}

}
\details{
The array is split into 8 sections and a $chi^2$ goodness of fit test is used to assess if the co-ordinates are equally spread among the 8 sections. The 
$chi^2$ statistic is defined such that a higher value indicates a less random distribution.

}

\value{
Numeric value giving the result of the $chi^2$ test.
}


\author{Mark Dunning}

\examples{

data(BLData)
co=getProbeCoords(BLData, 2, 1)
checkRandomness(BLData, 1, co)

}
\seealso{
\code{\link{getProbeCoords}}
}

\keyword{methods}
\name{plotDistances}
\alias{plotDistances}

\title{Plot the distances between a set of co-ordinates}

\description{
The function calculates the distance matrix for a set of co-ordinates and plots in the form of a histogram. The input co-ordinates could be the co-ordinates 
for a particular bead type on an array.
}

\usage{
plotDistances(coords)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{coords}{n by 2 vector where first column is list of x co-ordinates and second column is list of y co-ordinates
}
}
\details{

A histogram is produced with the mean of the distance matrix indicated in red. A normal distribution with mean = 650 and sd =350 is overlayed.

}

\value{
Plot is produced on the current graphical device
}


\author{Mark Dunning}

\examples{

data(BLData)
co=getProbeCoords(BLData, 2, 1)
plotDistances(co)

}
\seealso{
\code{\link{getProbeCoords}}

}
\keyword{methods}
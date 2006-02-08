\name{getProbeIDs}
\alias{getProbeIDs}

\title{Get unique ProbeIDs}
\description{
Return list of different ProbeIDs used on the array
}

\usage{
getProbeIDs(BLData)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{BLData}{\code{BeadLevelList} object containing bead-level data}
}
\value{
List of numeric values giving all ProbeIDs used on the array
}


\author{Mark Dunning}

\examples{

data(BLData)
getProbeIDs(BLData)

}

\keyword{methods}
\name{copyBeadLevelList}
\alias{copyBeadLevelList}
\title{Copies a BeadLevelList Object}
\description{
Make a new copy of a \code{BeadLevelList} object.
}
\usage{
copyBeadLevelList(object)
}
\arguments{
  \item{object}{\code{BeadLevelList}}
}
\details{
  \code{copyBeadLevelList} makes a new copy of a BeadLevelList object.
  This is necessary beacuse the beadData slot is stored as an
  environment.
  }
\value{
A new \code{\link[beadarray:BeadLevelList]{BeadLevelList}} object containing data from object.
}
\author{Matt Ritchie}
\examples{
data(BLData)
BLDataNew <- copyBeadLevelList(BLData)
BLData@beadData # the same bead level data is now
BLDataNew@beadData # stored in different environments
}
\keyword{manip}

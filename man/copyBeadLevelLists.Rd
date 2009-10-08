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
  \code{copyBeadLevelList} makes a new copy of a \code{BeadLevelList} object.
  This is necessary beacuse the beadData slot is stored as an environment.
  }
\value{
A new \code{\link[beadarray:class-BeadLevelList]{BeadLevelList}} object
containing the data from \code{object}.
}
\author{Matt Ritchie}
\examples{
data(BLData)
BLDataNew = copyBeadLevelList(BLData)
BLData@beadData # the same bead level data is now
BLDataNew@beadData # stored in different environments
}
\keyword{manip}

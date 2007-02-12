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
A \code{\link[beadarray:BeadLevelList]{BeadLevelList}} object holding data from all the arrays and all genes from the individual objects.
}
\author{Matt Ritchie}
\examples{
data(BLData)
BLDataNew <- copyBeadLevelList(BLData)
BLData@beadData # bead level data stored in different environments
BLDataNew@beadData
}
\keyword{manip}

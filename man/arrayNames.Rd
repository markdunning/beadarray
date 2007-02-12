\name{arrayNames}
\alias{arrayNames}
\title{Gets the array names from a BeadLevelList Object}
\description{
Retrieves the array names from a \code{BeadLevelList} object.
}
\usage{
arrayNames(object)
}
\arguments{
  \item{object}{\code{BeadLevelList}}
}
\details{
  \code{arrayNames} retrieves the name of each array from the arrayInfo
  slot.
   }
\value{
A character vector containing the names of the individual arrays.
}
\author{Matt Ritchie}
\examples{
data(BLData)
arrayNames(BLData)
}
\keyword{manip}

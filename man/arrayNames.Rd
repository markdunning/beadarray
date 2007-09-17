\name{arrayNames}
\alias{arrayNames}
\title{Gets the array names from a BeadLevelList Object}
\description{
Retrieves the array names from a \code{BeadLevelList} object.
}
\usage{
arrayNames(object, arrays=NULL)
}
\arguments{
  \item{object}{\code{BeadLevelList}}
  \item{arrays}{integer (scalar or vector) specifying the arrays to
    retrieve the names of.  When \code{NULL}, all the array names are returned.} 
}
\details{
  \code{arrayNames} retrieves the name of the array/s from the arrayInfo
  slot.
   }
\value{
A character vector containing the names of the individual array/s.
}
\author{Matt Ritchie}
\examples{
data(BLData)
arrayNames(BLData)
}
\keyword{manip}

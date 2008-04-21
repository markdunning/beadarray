\name{arrayNames}
\alias{arrayNames}
\title{Gets the strip/array names from a BeadLevelList Object}
\description{
Retrieves the strip/array names from a \code{BeadLevelList} object.
}
\usage{
arrayNames(object, arrays=NULL)
}
\arguments{
  \item{object}{\code{BeadLevelList}}
  \item{arrays}{integer (scalar or vector) specifying the strips/arrays to
    retrieve the names of.  When \code{NULL} the names of all strips/arrays
    are returned.} 
}
\details{
  \code{arrayNames} retrieves the name of the strip(s)/array(s) from the arrayInfo
  slot.
   }
\value{
A character vector containing the names of the individual strips(s)/array(s).
}
\author{Matt Ritchie}
\examples{
data(BLData)
arrayNames(BLData)
}
\keyword{manip}

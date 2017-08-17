\name{sectionNames}
\alias{sectionNames}
\title{Gets the section names from a beadLevelData Object}
\description{
Retrieves the section names from a \code{beadLevelData} object.
}
\usage{
sectionNames(object, arrays=NULL)
}
\arguments{
  \item{object}{Object of class \code{beadLevelData}}
  \item{arrays}{integer (scalar or vector) specifying the sections/arrays to retrieve the names of.  When \code{NULL} the names of all sections/arrays are returned.} 
}
\details{
  \code{sectionNames} retrieves the name of the sections from the sectionInfo slot.
   }
\value{
A character vector containing the names of the individual sections.
}
\author{Mark Dunning}
\examples{

if(require(beadarrayExampleData)){

data(exampleBLData)
sectionNames(exampleBLData)


}

}
\keyword{manip}

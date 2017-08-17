\name{getBeadData}
\alias{getBeadData}
\title{Get raw data from a beadLevelData object}
\description{
Retrieves the raw bead data from a \code{beadLevelData} object 
for a given section/array.
}
\usage{
getBeadData(BLData, what="Grn", array=1)
}
\arguments{
  \item{BLData}{\code{BeadLevelList}}
  \item{what}{character string specifying the values to retrieve (e.g. "ProbeID", "Grn" etc.).}
  \item{array}{integer specifying the section/array to use}
  }

\value{
A vector containing the specified bead data for the particular array.
}
\author{Mark Dunning}

\examples{

if(require(beadarrayExampleData)){

data(exampleBLData)

summary(getBeadData(exampleBLData))

}


}
\keyword{manip}

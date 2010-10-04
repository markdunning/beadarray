\name{getBeadData}
\alias{getBeadData}
\title{Get raw data from a beadLevelData object}
\description{
Retrieves the raw bead data from a \code{beadLevelData} object 
for a given strip/array.
}
\usage{
getBeadData(BLData, what="Grn", array=1)
}
\arguments{
  \item{BLData}{\code{BeadLevelList}}
  \item{what}{character string specifying the values to retrieve.}
  \item{array}{integer specifying the section/array to use}
  }

\value{
A vector containing the raw bead data (or residuals) for a particular array.
}
\author{Mark Dunning}
\examples{
data(BLData)
summary(getBeadData(BLData))
}
\keyword{manip}

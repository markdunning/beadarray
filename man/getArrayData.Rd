\name{getArrayData}
\alias{getArrayData}
\title{Gets the raw intensity data from a BeadLevelList object}
\description{
Retrieves the raw intensity information from a given array stored in a \code{BeadLevelList} object.
}
\usage{
getArrayData(BLData,which="G",array=1,log=TRUE)
}
\arguments{
  \item{BLData}{\code{BeadLevelList}}
  \item{which}{character string specifying which intensities to obtain.
    Possibilities are \code{"ProbeID"}, \code{"GrnX"}, \code{"GrnY"},
    \code{"G"}, \code{"Gb"} for single channel data and additional
    columns \code{"R"} and \code{"Rb"} for two-colour data}
  \item{array}{numeric value giving the array we are interested in}
  \item{log}{if TRUE log2 of the raw intensities are returned}
}
\details{
  \code{getArrayData} retrieves the raw data from \code{which} column of
  \code{array}.  When intensity data is chosen, the data can be on the
  log2 scale when \code{log==TRUE}.
}
\value{
A vector containing the intensities
}
\author{Matt Ritchie}
\examples{
data(BLData)
summary(getArrayData(BLData))
}
\keyword{manip}

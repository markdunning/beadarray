\name{getArrayData}
\alias{getArrayData}
\title{Gets the raw data from a BeadLevelList object}
\description{
Retrieves the per bead data from a \code{BeadLevelList} object 
for a given array.
}
\usage{
getArrayData(BLData,which="G",array=1,log=TRUE,n=3)
}
\arguments{
  \item{BLData}{\code{BeadLevelList}}
  \item{which}{character string specifying which intensities to obtain.
    Possibilities are \code{"ProbeID"}, \code{"GrnX"}, \code{"GrnY"},
    \code{"G"}, \code{"Gb"} for single channel data and \code{"R"},
    \code{"Rb"}, \code{"residR"}, \code{"residG"}, \code{"M"} (log-ratios)
    \code{"residM"} or \code{"A"} (average log-intensities) for
    two-colour data}
  \item{array}{numeric value giving the array we are interested in}
  \item{log}{if TRUE log2 of the raw intensities are returned}
  \item{n}{numeric value defining a cut-off of 'n' median absolute
    deviations from the median to determine outliers in
    \code{createBeadSummaryData}. The default value is 3.
    This argument is used to calculate the summarised R, G of M values
    which are used to calculate the residuals when \code{which="residR"}
    \code{"residG"} or \code{"residM"}.}
  }
\details{
  \code{getArrayData} retrieves the raw per bead data from a given
  array.  The data is either stored in \code{BeadLevelList} object
  (e.g. \code{"ProbeID"}, \code{"GrnX"}, \code{"GrnY"}, \code{"G"},
  \code{"Gb"}, \code{"R"} and \code{"Rb"} or calculated from these
  values (e.g. \code{"residR"}, \code{"residG"}, \code{"M"},
  \code{"residM"} or \code{"A"}).  When \code{log==TRUE}, intensity
  data is returned on the log2 scale.
}
\value{
A vector containing the bead data for a particular array.
}
\author{Matt Ritchie}
\examples{
data(BLData)
summary(getArrayData(BLData))
}
\keyword{manip}

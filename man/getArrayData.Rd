\name{getArrayData}
\alias{getArrayData}
\title{Gets the raw data from a BeadLevelList object}
\description{
Retrieves the per bead data from a \code{BeadLevelList} object 
for a given array.
}
\usage{
getArrayData(BLData, what="G", array=1,log=TRUE, method="illumina", n=3, trim=0.05)
}
\arguments{
  \item{BLData}{\code{BeadLevelList}}
  \item{what}{character string specifying the values to retrieve.
    Possibilities are \code{"ProbeID"}, \code{"GrnX"}, \code{"GrnY"},
    \code{"G"}, \code{"Gb"} for single channel data and \code{"R"},
    \code{"Rb"}, \code{"residR"}, \code{"residG"}, \code{"M"} (log-ratios)
    \code{"residM"}, \code{"A"} (average log-intensities) and \code{"beta"} 
    (=R/(R+G)) for two-colour data}
  \item{array}{numeric value giving the array we are interested in}
  \item{log}{if TRUE log2 of the raw intensities are returned (ignored if \code{what="beta"})}
  \item{method}{chracter string specifying the summarisation method to
    use in \code{createBeadSummaryData} (see help page for further
    details). Only used when \code{what="residR"}, \code{"residG"}
    or \code{"residM"}.}
  \item{n}{use n median absolute deviations from the median as a
    cut-off for outliers. Only used when \code{what="residR"}, \code{"residG"}
    or \code{"residM"} and \code{method="illumina"}.}
  \item{trim}{fraction of intensities to remove from the bead summary
    calculations when \code{method="trim"}, or the fraction of
    intensities to set to the \code{trim} and 1-\code{trim} percentile
    intensities when \code{method="winsorize"}.  Default value is 0.05.
    Only used when \code{what="residR"}, \code{"residG"} or \code{"residM"}.}
  }
\details{
  \code{getArrayData} retrieves the raw per bead data from a given
  array.  The data is either extracted from a \code{BeadLevelList} object
  (e.g. \code{"ProbeID"}, \code{"GrnX"}, \code{"GrnY"}, \code{"G"},
  \code{"Gb"}, \code{"R"} and \code{"Rb"} or calculated from these
  values (e.g. \code{"residR"}, \code{"residG"}, \code{"M"},
  \code{"residM"}, \code{"A"} or \code{"beta"}).
  When \code{log=TRUE}, intensity data is returned on the log2 scale.
}
\value{
A vector containing the raw bead data (or residuals) for a particular array.
}
\author{Matt Ritchie}
\examples{
data(BLData)
summary(getArrayData(BLData))
}
\keyword{manip}

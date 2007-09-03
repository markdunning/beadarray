\name{beadResids}
\alias{beadResids}
\title{Calculates per array bead-level residuals}
\description{
Retrieves the per bead data stored in a \code{BeadLevelList} object
for a given array.  In the case of the log-ratios or average
log-intensities, these values are calculated for each bead.
}
\usage{
beadResids(BLData, what="G", array=1, log=TRUE,
                 method="illumina", n=3, trim=0.05)
}
\arguments{
  \item{BLData}{\code{BeadLevelList}}
  \item{what}{character string specifying which intensities to obtain.
    Possibilities are \code{"G"} (single channel data) and \code{"R"}
    (two-colour data)}
  \item{array}{numeric value giving the array we are interested in}
  \item{log}{if TRUE then use log2 intensities of each bead}
  \item{method}{chracter string specifying the summarisation method 
    (see help page for \code{createBeadSummaryData} for further details).}
  \item{n}{use n median absolute deviations from the median as a
    cut-off for outliers. Only used when \code{what="residR"}, \code{"residG"}
    or \code{"residM"} and \code{method="illumina"} (see
    \code{createBeadSummaryData} help page for further details).}
  \item{trim}{fraction of intensities to remove from the bead summary
    calculations when \code{method="trim"}, or the fraction of
    intensities to set to the \code{trim} and 1-\code{trim} percentile
    intensities when \code{method="winsorize"}.  Default value is 0.05.
    Only used when \code{what="residR"}, \code{"residG"} or \code{"residM"}.}
  }
\details{
  Calculates the residuals, i.e. the differences between the summary
  values obtained from \code{createBeadSummaryData} and the individual
  values for each bead.
}
\value{
A vector containing the residual values.
}
\author{Matt Ritchie}
\examples{
data(BLData)
summary(beadResids(BLData, log=TRUE))
}
\keyword{manip}

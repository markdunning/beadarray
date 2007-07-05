\name{beadResids}
\alias{beadResids}
\title{Calculates per array bead-level residuals}
\description{
Retrieves the per bead data stored in a \code{BeadLevelList} object
for a given array.  In the case of the log-ratios or average
log-intensities, these values are calculated for each bead.
}
\usage{
beadResids(BLData,which="G",array=1,log=TRUE,n=3)
}
\arguments{
  \item{BLData}{\code{BeadLevelList}}
  \item{which}{character string specifying which intensities to obtain.
    Possibilities are \code{"G"} (single channel data) and \code{"R"}
    (two-colour data)}
  \item{array}{numeric value giving the array we are interested in}
  \item{log}{if TRUE log2 of the intensities are used}
  \item{n}{numeric value defining a cut-off of 'n' median absolute
    deviations from the median to determine outliers in
    \code{createBeadSummaryData}. The default value is 3.
    This argument is used to calculate the summarised R or G values.}
  }
\details{
  Calculates the residuals, i.e. the differences between the summary
  values obtained from \code{createBeadSummaryData} and the individual
  values for each bead.  The residuals can be calculated on the original
  or log2 scale.
}
\value{
A vector containing the residual values.
}
\author{Matt Ritchie}
\examples{
data(BLData)
summary(beadResids(BLData))
}
\keyword{manip}

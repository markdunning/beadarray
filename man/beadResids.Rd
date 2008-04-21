\name{beadResids}
\alias{beadResids}
\title{Calculates per strip/array bead-level residuals}
\description{
Calculates the per bead residuals for a given strip/array using 
data from a \code{BeadLevelList}.
}
\usage{
beadResids(BLData, what="G", array=1, log=TRUE,
                 method="illumina", n=3, trim=0.05)
}
\arguments{
  \item{BLData}{\code{BeadLevelList}}
  \item{what}{character string specifying which intensities to use in the 
    calculation of residuals.  See \code{getArrayData} for a list of possibilities}
  \item{array}{integer specifying the strip/array to use}
  \item{log}{if TRUE then use log2 intensities of each bead}
  \item{method}{character string specifying the summarisation method 
    (see help page for \code{createBeadSummaryData} for further details).}
  \item{n}{numeric value defining a cut-off for the number of 
    median absolute deviations (MADs) from the median to use for 
    determining outliers.  The default value is 3.  Only used when 
    \code{method="illumina"} (see \code{createBeadSummaryData} help 
    page for further details).}
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

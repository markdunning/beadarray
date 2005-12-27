\name{probeDiagnostics}
\alias{probeDiagnostics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Diagnostics plot for a bead type}
\description{
Function to produce 4 diagnostic plots on the same page for a particular bead type on an array
}
\usage{
probeDiagnostics(BLData, probe, array, label = TRUE,
log = FALSE, n = 3, showOutliers = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{BeadLevelList object}
  \item{probe}{numeric value giving the probeID of the bead type we want to plot}
  \item{array}{numeric value of the array we want}
  \item{label}{if TRUE all beads will be assigned a numeric label}
  \item{log}{if TRUE outliers for this bead type will be calculated using log2 bead intensities}
  \item{n}{outliers for this bead type are beads with intensities more
    than 'n' standard deviations from the mean}
  \item{showOutliers}{boolean specifying whether the outliers should be
    shown on the plots}
\item{\dots}{further graphical parameters to \code{plot}}
}
\details{
Combines the results of 4 seperate plotting functions onto the same page
}
\value{
Plot is produced on the current graphical device.

}

\author{Mark Dunning}


\examples{

data(BLData)

probeDiagnostics(BLData, 2, 1)

probeDiagnostics(BLData, 2, 1, log=TRUE)

probeDiagnostics(BLData, 2, 1, log=TRUE, label=FALSE)

}
\keyword{hplot}
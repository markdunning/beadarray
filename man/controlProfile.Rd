\name{controlProfile}
\alias{controlProfile}
\docType{data}
\title{Control annotation for Illumina expression chips}
\description{Data frame defining control information for the example \code{\link{BLData}} object included with the beadarray package. The first column contains ArrayAddressIDs to be found in the \code{beadLevelData} object and which control type they represent.}
\usage{data(controlProfile)}

\examples{
library(beadarray)
data(controlProfile)
head(controlProfile)
table(controlProfile)

}

\keyword{datasets}

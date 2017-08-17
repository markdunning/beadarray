\name{medianNormalise}
\alias{medianNormalise}
\title{Median normalise data in a matrix}
\description{
  Normalises expression intensities so that the intensities or
  log-ratios have equal median values across a series of arrays (columns).
}
\usage{
medianNormalise(exprs, log=TRUE)
}
\arguments{
  \item{exprs}{a matrix of expression values}
  \item{log}{if TRUE then do a log2 transformation prior to normalising}
}
\details{
  Normalisation is intended to remove from the expression measures any
  systematic trends which arise from the microarray technology rather
  than from differences between the probes or between the target RNA
  samples hybridized to the arrays.

  For median normalisation, the intensity for each gene is adjusted by
  subtracting the median of all genes on the array and then adding the
  median across all arrays. The effect is that each array then has the
  same median value.
}
\value{
  Produces a matrix of normalised intensity values (on the log2 scale by
  default) with the same dimensions as \code{exprs}.
}
\author{Mark Dunning}

\examples{

if(require(beadarrayExampleData)){

data(exampleSummaryData)

exampleSummaryData.log2 <- channel(exampleSummaryData,"G")

exampleSummaryData.med = assayDataElementReplace(exampleSummaryData.log2, "exprs", medianNormalise(exprs(exampleSummaryData.log2)))

}

}
\keyword{methods}
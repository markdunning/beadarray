\name{medianNormalise}
\alias{medianNormalise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Median Normalisation for BeadLevelList objects }
\description{
Normalizes expression intensities so that the intensities or
     log-ratios have equal median values across a series of arrays.
}
\usage{
medianNormalise(BLData)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{can be either an BeadLevelList object or an BeadSummaryList object.} 
}
\details{
Normalization is intended to remove from the expression measures any systematic trends which arise from the microarray technology rather than from differences between the probes or between the target RNA samples hybridized to the arrays.

For median normalisation, the intensity for each gene is adjusted by subtracting the median of all genes on the array and then adding the median across all arrays. The effect is that each array then has the same median value.

}

\value{
Produces an BeadLevelList containing normalised intensity values. Values will be on the 
log2 scale

}
\author{Mark Dunning}


\examples{

data(BLData)

BLData.med = medianNormalise(BLData)

#All distributions are fixed to have same median

boxplot(BLData.med$G~col(BLData.med$G))

}

\keyword{methods}
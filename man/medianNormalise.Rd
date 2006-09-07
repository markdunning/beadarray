\name{medianNormalise}
\alias{medianNormalise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Median Normalisation for BeadLevelList objects }
\description{
Normalizes expression intensities so that the intensities or
     log-ratios have equal median values across a series of arrays.
}
\usage{
medianNormalise(exprs, log=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exprs}{can be either the G matrix from a BeadLevelList object or matrix of expression values} 
  \item{log}{if TRUE then do a log2 transformation prior to normalising}
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

data(BSData)

BSData.med = assayDataElementReplace(BSData, "exprs", medianNormalise(exprs(BSData)))




}

\keyword{methods}
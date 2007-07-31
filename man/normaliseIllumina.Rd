\name{normaliseIllumina}
\alias{normaliseIllumina}
\title{Normalise Illumina expression data}
\description{
Normalises expression intensities so that the intensities are comparable between arrays.
}
\usage{
normaliseIllumina(BSData, method="quantile", transform="none", T=NULL)
}

\arguments{
  \item{BSData}{an 'ExpressionSetIllumina' object}
  \item{method}{character string specifying normalisation method
    (options are \code{"quantile"}, \code{"qspline"}, \code{"vsn"},
    \code{"rankInvariant"}, \code{"median"} and \code{"none"}.}
  \item{transform}{character string specifying transformation to apply
    to the data (options are \code{"none"}, \code{"log2"} and \code{"vst"}}
  \item{T}{A target distribution vector used in \code{"rankInvariant"}
    normalisation. If \code{NULL}, the mean is used.}  
}
\details{
  Normalisation is intended to remove from the expression measures any
  systematic trends which arise from the microarray technology rather
  than from differences between the probes or between the target RNA
  samples hybridized to the arrays.


  When \code{transform="vst"} the variance stabilising transformation
  from the 'lumi' package is applied to the data.  Refer to the 'lumi'
  documentation for further particulars.

  For further particulars on the different normalisation methods options
  refer to the individual help pages (\code{?normalize.quantiles} for
  \code{"quantile"}, \code{?normalize.qspline} for "qspline",
  \code{?rankInvariantNormalise} for \code{"rankInvariant"},
  \code{?medianNormalise} for \code{"median"} and \code{?vsn2}
  for \code{"vsn"}.

  For median normalisation, the intensity for each gene is adjusted by
  subtracting the median of all genes on the array and then adding the
  median across all arrays. The effect is that each array then has the
  same median value.

  Note: If your \code{BSData} object contains data already on the
  log-scale, be careful that you choose an appropriate \code{transform}
  to avoid logging it twice.  The same applies for the \code{"vst"}
  transformation and \code{"vsn"} normalisation methods which require
  the expression data stored in \code{BSData} to be on the original
  (un-logged) scale.
}

\value{
  An 'ExpressionSetIllumina' object which conatains the transformed
  and normalised expression values for each array.
}

\author{Matt Ritchie}

\examples{
data(BSData)
BSData.norm = normaliseIllumina(BSData, method="quantile", transform="log2")
}

\keyword{methods}
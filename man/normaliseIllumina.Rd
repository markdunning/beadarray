\name{normaliseIllumina}
\alias{normaliseIllumina}
\title{Normalise Illumina expression data}
\description{
Normalises expression intensities from an \code{ExpressionSetIllumina} 
object so that the intensities are comparable between arrays.
}
\usage{
normaliseIllumina(BSData, method="quantile", transform="none", T=NULL, status=fData(BSData)$Status,negctrl="negative",regular="regular",...)
}

\arguments{
  \item{BSData}{an \code{ExpressionSetIllumina} object}
  \item{method}{character string specifying normalisation method
    (options are \code{"quantile"}, \code{"qspline"}, \code{"vsn"},
    \code{"rankInvariant"}, \code{"median"} and \code{"none"}.}
  \item{transform}{character string specifying transformation to apply
    to the data prior to normalisation (options are \code{"none"}, 
    \code{"log2"}, \code{neqc}, \code{rsn} and \code{"vst"}}
  \item{T}{A target distribution vector used when \code{method="rankInvariant"}
    normalisation.  If \code{NULL}, the mean is used.}

  \item{status}{character vector giving probe types (used in neqc normalisation only)		
}

\item{negctrl}{
character vector giving negative control probes (used in neqc normalisation only)
}

\item{regular}{

character vector giving regular probes (used in neqc normalisation only)
}


  \item{...}{further arguments to be passed to \code{lumiT} or \code{neqc}}
}
\details{
  Normalisation is intended to remove from the expression measures any
  systematic trends which arise from the microarray technology rather
  than from differences between the probes or between the target RNA
  samples hybridized to the arrays.
 
  In this function, the \code{transform} specified by the user is 
  applied prior to the chosen normalisation procedure.

  When \code{transform="vst"} the variance-stabilising transformation
  from the 'lumi' package is applied to the data.
  Refer to the \code{lumiT} documentation for further particulars.  Note that 
  the Detection P values are only passed on when they are available 
  (i.e. not NA). The rsn option calls code directly from \code{lumi}.

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
  to avoid transforming it twice.  The same applies for the \code{"vst"}
  transformation and \code{"vsn"} normalisation methods which require
  the expression data stored in \code{BSData} to be on the original
  (un-logged) scale.  When \code{method="vsn"}, transform must be set to 
  \code{"none"}, since this method transforms and normalises the data as 
  part of the model.


  The neqc normalisation is described in Shi et al (2010) and documented in the limma package. Note that the output from this method has control probes removed.
	
}

\value{
  An 'ExpressionSetIllumina' object which conatains the transformed
  and normalised expression values for each array.
}

\author{Matt Ritchie and Mark Dunning}

\examples{

if(require(beadarrayExampleData)){

data(exampleSummaryData)
exampleSummaryData.norm = normaliseIllumina(channel(exampleSummaryData, "G"), method="quantile", transform="none")

exampleSummaryData.rsn = normaliseIllumina(channel(exampleSummaryData, "G.ul"),method="rsn", transform="none")

}

}

\references{
Shi. W, Oshlack. A, Smyth, GK, (2010) Optimizing the noise versus bias trade-off for Illumina whole genome expression BeadChips. Nucleic Acids Research

Lin, S.M., Du, P., Kibbe, W.A., (2008) 'Model-based Variance-stabilizing Transformation for Illumina Microarray Data', Nucleic Acids Res. 36, e11

}

\keyword{methods}

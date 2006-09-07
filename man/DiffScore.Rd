\name{DiffScore}
\alias{DiffScore}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Illumina's DE statistic}
\description{
Function to compute the DiffScore statistic to assess the probability of differential expression of genes between a reference and condition array.

NB At present this function does not compute the same values as output by BeadStudio


}
\usage{
DiffScore(BSData, QC, cond, ref)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BSData}{ExpressionSetIllumina object representing bead summary data}
  \item{QC}{quality control information created by readQC}
  \item{cond}{index for the condition array}
  \item{ref}{index for the reference array}	

}

\details{
The model assumes that bead variation is comprised or nonspecific biological variation and technical error. Technical error
is estimated using a robust least squares fit of the standard error and mean for beads on the reference and condition arrays respectively. Biological 
variation is estimated using the average intensity of the negative controls on the reference and condition arrays (found in the QC object). A p-value for
the statistic is computed by dividing the absolute difference in expression on reference and condition arrays by the sum of the technical and biological varation on 
condition and reference arrays. The DiffScore is then

$10*sign(mu_cond - mu_ref)log_10(p)$

Note that this statistic is intended to be used on unlogged intensities. For log$_2$ data, the lmFit function from limma can be easily used on the exprs matrix.



For more details see the BeadStudio manual and "A model of Technical Variation of Microarray Signal" - Chudin et al. Journal of Comp Bio. Vol 13, 4 (2006)


}
\value{
vector of DiffScore values, one for each gene in the BSData object.
}
\seealso{\code{\link[limma:lmFit]{lmFit}}}
\author{Mark Dunning}

\keyword{IO}

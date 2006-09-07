\name{backgrondNormalise}
\alias{backgroundNormalise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Background Normalisation for Bead Summary Data objects }
\description{
Normalise expression matrix using values of negative controls
}
\usage{
backgroundNormalise(exprs, QC)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exprs}{An expression matrix from a ExpressionSet object} 
  \item{QC}{A quality control assayData object created by readQC}
}
\details{

This function implements the background normalisation method recommened by Illumina. All normalisation options in BeadStudio (with the exception of no normalisation) include background normalisation as a first step.

The QC object contains averaged values for the Illumina negative controls placed on each array in an experiment. For background normalisation, we simply 
subtract the value of the negative control on a particular array from all bead summary values on that array.

Note that this background normalisation is in addition to a local background correction performed at the bead level and is liable to create negative values for genes with low expression. 
It is therefore not advisable if you wish to analyse data on a log2 scale.

}


\value{
Expression matrix containing normalise values with the same dimensions as the input

}
\author{Mark Dunning}


\examples{

data(BSData)
data(QC)

BSData.bgnorm = assayDataElementReplace(BSData, "exprs", backgroundNormalise(exprs(BSData), QC))



}

\seealso{\code{\link{readQC}}}

\keyword{methods}
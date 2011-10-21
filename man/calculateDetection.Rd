\name{calculateDetection}
\alias{calculateDetection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate detection scores}
\description{
Function to calculate detection scores for summarized data if they are not available.
}
\usage{
calculateDetection(BSData, status=fData(BSData)$Status, negativeLabel="negative")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BSData}{An ExpressionSetIllumina object}
  \item{status}{character vector giving probe types}
  \item{negativeLabel}{character giving identifer for negative controls}
 
}
\details{

Detection scores are a measure of whether the probe is showing any specific expression. This function implements Illumina's method for calculating the detection scores for all bead types on a given array. Within an array, Illumina discard negative control bead-types whose summary values are more than three MADs from the median for the negative controls. Illumina then rank the summarized intensity for each other bead-type against the summarized values for the remaining negative control bead-types and calculate a detection p-value  1-R/N, where R is the relative rank of the bead intensity when compared to the $N$ remaining negative controls. Thus, if a particular bead has higher intensity than all the negative controls it will be assigned a value of 0. This calculation is repeated for all arrays. \\

The function expects the negative controls to be indicated by the \code{Status} column in the \code{featureData} slot of the \code{ExpressionSetIllumina} object. If this is not present the user can supply a status vector with the same length as the number of rows in the ExpressionSetIllumina object.
}
\value{
Matrix of detection scores with the same dimensions as the exprs matrix of BSData. This matrix can be stored in a ExpressionSetIllumina object using the \code{Detection} function
}

\examples{

if(require(beadarrayExampleData)){

data(exampleSummaryData)

table(fData(exampleSummaryData)[,"Status"])

exampleSummaryData.log2 <- channel(exampleSummaryData ,"G")

det <- calculateDetection(exampleSummaryData.log2)

Detection(exampleSummaryData.log2) <- det

}

}

\author{Mark Dunning and Andy Lynch}

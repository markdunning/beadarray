\name{calculateDetection}
\alias{calculateDetection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate detection scores }
\description{
Function to calculate detection scores for summarized data if they are not available.
}
\usage{
calculateDetection(BSData)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BSData}{An ExpressionSetIllumina object}
}
\details{

The function implements Illumina's method for calculating the detection scores for all bead types on a given array. 
Within an array, Illumina discard negative control bead-types whose summary values are more than three MADs from the median for the negative controls. Illumina then rank the summarized intensity for each other bead-type against the summarized values for the remaining negative control bead-types and calculate a detection p-value  1-R/N, where R is the relative rank of the bead intensity when compared to the $N$ remaining negative controls. Thus, if a particular bead has higher intensity than all the negative controls it will be assigned a value of 0. 
This calculation is repeated for all arrays stored in the BSData object.
}
\value{
A modified version of BSData is returned with the detection slot filled with the newly created scores.
}
\author{Mark Dunning and Andy Lynch}

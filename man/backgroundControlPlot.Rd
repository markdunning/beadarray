\name{backgroundControlPlot}
\alias{backgroundControlPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{QA measures based on bead-level negative controls}
\description{
Function for plotting the bead-level intensities for all the negative controls that are placed on an array. Typically there are around a thousand of these controls, each replicated 30 times. 
The sequences used for these controls should not target any part of the genome and therefore we should not observe any signal. 

}
\usage{
backgroundControlPlot(BLData, array = 1, plot = FALSE, t1,t2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{A BeadLevelList for an Illumina expression chip}
  \item{array}{The number of the array in BLData that we want QA of.}
  \item{plot}{if TRUE a diagnostic plot will be produced}
  \item{t1}{Pre-computed list of array intensities}
  \item{t2}{List of ProbeIDs for all beads on the array}
}
\details{
For QA, we report the mean and variance of all negative controls (of all bead-types) after first removing outiers using a 3 MAD cut-off. To retrieve the IDs of the negative controls, we make use of the annotation slot stored with the BeadLevelList object. It is therefore important that this information is accurate. 
A plot of all negative control bead-types can also be produced, where each bead-type is represented by a vertical line covering the inter-quartile range and ordered according to mean intensity. Too many high intensity values for the negatives could indicate a poor quality array.
}
\value{
The function returns the mean (AveNeg) and variance (VarNeg) of all negative control beads and a diagnositc plot if requested.
}
\author{Mark Dunning and Andy Lynch}
\seealso{\code{\link{calculateBeadLevelScores}}}


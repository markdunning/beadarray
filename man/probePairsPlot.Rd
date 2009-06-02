\name{probePairsPlot}
\alias{probePairsPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{QA plot using perfect match and mismatch controls}
\description{
QA using the hybridisations controls that are identical expect for a few bases. Hence lower signal should be observed for the mismatch probes.
}
\usage{
probePairsPlot(BLData, array = 1,t1,t2,plot=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{A BeadLevelList object for an Illumina expression array}
  \item{array}{Number of the array in BLData that we want QA for}
  \item{t1}{Pre-computed list of array intensities}
  \item{t2}{List of ProbeIDs for all beads on the array}
  \item{plot}{Return the plot or not}
}
\details{

Firstly, the annotation of the control probes found on the BeadLevelList object is retrieved. We then search the annoation for the perfect and mismatch controls (labeled phage\_lambda\_genome:pm and phage\_lambda\_genome:mm2). 
For each perfect match we find the corresponding mismatch control by comparing the sequences and plot the perfect match control next to the mismatch control.

}


\value{
Plot on the current graphical device with the perfect match beads in red and mismatches in purple. The labels on the x-axis show the ID of the controls. 
}
\author{Mark Dunning}
\seealso{\code{\link{calculateBeadLevelScores}}}

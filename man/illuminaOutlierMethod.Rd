\name{illuminaOutlierMethod}
\alias{illuminaOutlierMethod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identifier outliers on an array section
}
\description{
Implementation of the illumina method for excluding outliers using a fixed number of MADs (median absolute deviations) cutoff for each bead type
}
\usage{
illuminaOutlierMethod(inten, probeList, wts=1, n = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inten}{
a list of intensities
}
  \item{probeList}{
the IDs corresponding to each intensity value
}
  \item{wts}{
Weights associated with beads, indicating those recommended for removal by, for example, \code{\link{BASH}}
}

\item{n}{
number of MADs cutoff used
}
}

\details{
This function is called within the summarisation routine of beadarray to exclude outliers from an array-section prior to summary. The intensities are not assumed to be on any particular scale and can result from any user-defined transformation function. 

Beads with weight zero do not contribute to the outlier calling.
}


\value{
the positions in the original vector that were determined to be outliers

}
\author{
Mark Dunning
}

\seealso{\code{\link{squeezedVarOutlierMethod}}}


\examples{

if(require(beadarrayExampleData)){

data(exampleBLData)

oList = illuminaOutlierMethod(logGreenChannelTransform(exampleBLData, 1), getBeadData(exampleBLData, array=1, what="ProbeID"))
 
}

}


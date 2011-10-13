\name{weightsOutlierMethod}
\alias{weightsOutlierMethod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
returns the locations of all beads with weight=0.
}
\description{
This function identifies those beads that have been set to having a weight of 0. Primarily intended to allow deprecation of the function showArrayMask.Rd by the more flexible outlierplot.
}
\usage{
weightsOutlierMethod(inten, probeList, wts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inten}{
a list of intensities but the intensities are not in fact used
}
  \item{probeList}{
the IDs corresponding to each intensity value (not used)
}
  \item{wts}{
Weights associated with beads. 
}
}

\details{
This function the locations at which wts is equal to zero.
}


\value{
the positions in the original vector that had weight zero

}
\author{
Andy Lynch
}

\seealso{\code{\link{squeezedVarOutlierMethod}}}


\examples{


if(require(beadarrayExampleData)){
data(exampleBLdata)
outlierplot(exampleBLdata, array=1, outlierFun = weightsOutlierMethod,horizontal = FALSE)
}


}


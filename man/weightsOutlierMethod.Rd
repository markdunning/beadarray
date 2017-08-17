\name{weightsOutlierMethod}
\alias{weightsOutlierMethod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
returns all beads with weight=0.
}
\description{
This function identifies those beads that have been set to having a weight of 0. Primarily intended to allow deprecation of the function \code{showArrayMask} by the more flexible outlierplot.
}
\usage{
weightsOutlierMethod(inten, probeList, wts,n=3)
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
  \item{n}{
another parameter that is not, in fact, used. 
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

\dontrun{

if(require(beadarrayExampleData)){
data(exampleBLData)
outlierplot(exampleBLData, array=1, outlierFun = weightsOutlierMethod,horizontal = FALSE)
}

}

}


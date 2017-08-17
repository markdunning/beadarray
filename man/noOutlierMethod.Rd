\name{noOutlierMethod}
\alias{noOutlierMethod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
returns no outliers on an array section
}
\description{
The \code{summarize} function demands that an outlier function is called. This function allows one to satisfy this requirement and still not remove any outliers.
}
\usage{
noOutlierMethod(inten, probeList, wts=1,n=3)
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
Weights associated with beads. Again, these are not actually used.
}
  \item{n}{
another parameter that is not, in fact, used. 
}
}

\details{
This function returns integer(0).
}


\value{
integer(0)

}
\author{
Andy Lynch
}

\seealso{\code{\link{squeezedVarOutlierMethod}}}


\examples{

if(require(beadarrayExampleData)){

data(exampleBLData)

oList = noOutlierMethod(logGreenChannelTransform(exampleBLData, 1), getBeadData(exampleBLData, array=1, what="ProbeID"))
 
}

}


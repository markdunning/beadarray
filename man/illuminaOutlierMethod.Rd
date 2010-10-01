\name{illuminaOutlierMethod}
\alias{illuminaOutlierMethod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identifier outliers on an array section
}
\description{
Implementation of the illumina method for excluding outliers using a 3 MAD (median absolute deviation) cutoff for each bead type
}
\usage{
illuminaOutlierMethod(inten, probeList, n = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inten}{
a list of intensities
}
  \item{probeList}{
the IDs corresponding to each intensity value
}
  \item{n}{
number of MADs cutoff used
}
}
\value{
the positions in the original vector that were determined to be outliers

}
\author{
Mark Dunning
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(BLData)

oList = illuminaOutlierMethod(logGreenChannelTransform(BLData, 1), getBeadData(BLData, array=1, what="ProbeID"))
 

}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

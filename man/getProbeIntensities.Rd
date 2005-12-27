\name{getProbeIntensities}
\alias{getProbeIntensities}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get Intensites of Beads}
\description{
 Function to return the intensity of beads on an array with a certain probeID.
}
\usage{
getProbeIntensities(BLData, probe, array, log = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{BeadLevelList object containing bead level data}
  \item{probe}{numeric value for the probeID we want values for}
  \item{array}{numeric value giving the array we want to find values on}
  \item{log}{if TRUE then log2-transformed data will be returned}
}
\details{
  The function simply looks up the BLData\$R values for beads whose BLData\$probeID matches the given 'probeID'. A log2 transform is then applied if 
specified.
}

\value{
  vector of numeric values corresponding to the intensity of all beads on the array with the specified 'probeID'
}

\author{Mark Dunning}

\examples{

data(BLData)

getProbeIntensities(BLData, 2, 1)

##To return log2 data use

getProbeIntensities(BLData, 2, 1, log=FALSE)

#ProbeID 4279 is a housekeeping control in our experiment

getProbeIntensities(BLData, 4279, 1)

}
\keyword{methods}

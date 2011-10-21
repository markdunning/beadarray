\name{identifyControlBeads}
\alias{beadStatusVector}
\alias{identifyControlBeads}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Classify each bead according to its control status
}
\description{
Using the control annotation specified for the array, the function will classify each bead as belonging to a control group, or as being a regular probe.
}
\usage{
identifyControlBeads(BLData, array = 1, controlProfile = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{
a \code{beadLevelData} object
}
  \item{array}{
the numeric id of the array section
}
  \item{controlProfile}{
an optional control profile data frame
}
}
\details{
The function requires either that a control profile data frame is specified (This associates probe IDs with their control status - see the example), or that the annotation of the \code{beadLevelData} object be set to an array that the package recognises (see \code{\link{getAnnotation}}). Note that some positive control bead-types may also be functioning 'regular' probes.
}
\value{
a vector of character strings giving the status of each bead
}
\author{
Mark Dunning
}


\examples{

if(require(beadarrayExampleData)){

data(exampleBLData)

statVec = identifyControlBeads(exampleBLData)

table(statVec)


}

}


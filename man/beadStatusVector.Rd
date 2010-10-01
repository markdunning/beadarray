\name{beadStatusVector}
\alias{beadStatusVector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Classify each bead according to its control status
}
\description{
Using the control annotation specified for the array, the function will classify each bead as belonging to a control group, or being a regular probe.
}
\usage{
beadStatusVector(BLData, array = 1, controlProfile = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{
a \code{beadLevelData} object
}
  \item{array}{
the numeric id of the section
}
  \item{controlProfile}{
an optional control profile data frame
}
}
\details{
The function requires a control profile data frame which must be specified or created using the annotation of the \code{beadLevelData} object. 
}
\value{
a vector of character strings giving the status of each bead

}
\author{
Mark Dunning
}


\examples{

data(BLData)

data(controlProfile)

statVec = beadStatusVector(BLData, controlProfile=controlProfile)

table(statVec)

}


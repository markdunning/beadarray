\name{quickSummary}
\alias{quickSummary}
\title{
Create summary values for specified IDs
}
\description{
A utility function for quickly creating summary values for particular IDs (e.g. control IDs) on a given section.
}
\usage{
quickSummary(BLData, array = 1, transFun = logGreenChannelTransform, reporterIDs = NULL, reporterTags = NULL, reporterFun = function(x) mean(x, na.rm = TRUE))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{
a \code{beadLevelData} object
}
  \item{array}{
Which section to summarize
}
  \item{transFun}{
a transformation to be applied prior to summarization
}
  \item{reporterIDs}{
vector specifying the set of IDs to be summarized
}
  \item{reporterTags}{
a vector that divides the supplied IDs into categories
}
  \item{reporterFun}{
a function used to summarize each category
}
}
\details{
The function can be used to calculate summarized values for particular control types on a section. The IDs for all controls are supplied in the \code{reporterIDs} argument along with which control type they belong to in the \code{reporterTags} argument. A summarized value for each control type is then calculated with the specified function (default is mean).

}

\examples{

if(require(beadarrayExampleData)){

data(exampleBLData)

quickSummary(exampleBLData, array=1)
quickSummary(exampleBLData, array=2)

}

}

\author{
Mark Dunning
}



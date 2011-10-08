\name{addFeatureData}
\alias{addFeatureData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add probe data
}
\description{

Adds extra probe-specific data to an ExpressionSetIllumina object. 

}
\usage{
addFeatureData(data, toAdd = c("SYMBOL", "PROBEQUALITY", "CODINGZONE", "PROBESEQUENCE"), annotation = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
An ExpressionSetIllumina object
}
  \item{toAdd}{
Either a pre-prepared data frame, or characters which refer to mappings within an annotation pacakge
}
  \item{annotation}{
Optional character identifying the annotation of the ExpressionSetIllumina object. e.g. Humanv3, Mousev2. 
}
}
\details{

The function will identify which package should be used by concatenating illumina with the value of the annotation slot of the object, or the annotation argument passed in. If this package is not installed on the users computer, then the function will fail.

Assuming the package has been correctly loaded, the character vector toAdd is converted to the names of environments within the package. These environments are then queried with the featureNames of the input object. The result of each query is converted to a data frame and merged with the original feature data of the object.

Alternatively, rather than querying from an annotation pacakge, a pre-prepared data frame can be used.


}

\value{
An ExpressionSetIllumina object with modified featureData
}
\examples{

if(require(beadarrayExampleData)){

data(exampleSummaryData)

exampleSummaryData <- addFeatureData(exampleSummaryData)

}

}
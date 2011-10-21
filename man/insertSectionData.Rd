\name{insertSectionData}
\alias{insertSectionData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Modify the sectionData slot
}
\description{
A function to modify the sectionData slot of a \code{beadLevelData} object. Data can be be added if it is a data frame with a number of rows equal to the number of sections in the beadLevelData object.

}
\usage{
insertSectionData(BLData, what, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{
a \code{beadLevelData} object
}
  \item{what}{
a character string specifiying a name for the new data
}
  \item{data}{
a data frame containing the data we wish to add
}
}
\details{
This function allows users to modifiy the per\-section information that is included in the sectionData slot. Typical usage would be to store quality control data that has been computed.
}
\value{
a modified \code{beadLevelData} object with the new data attached to sectionData

}
\author{
Mark Dunning
}
\examples{

if(require(beadarrayExampleData)){

data(exampleBLData)

qct = makeQCTable(exampleBLData)

exampleBLData = insertSectionData(exampleBLData, what="ProbeQC", data = qct)

exampleBLData@sectionData

}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.


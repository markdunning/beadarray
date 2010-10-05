\name{ExpressionControlData}
\alias{ExpressionControlData}
\docType{data}
\title{Control annotation for Illumina expression chips}
\description{Data frames derived from the bgx files from Illumina that give details of the control probes used on Illumina expression arrays. A list structure is used with the control probes for a particular platform accessed by name. Note that the HumanHT12 arrays use the same probes and the Humanv3 and therefore the same annotation can be used.}
\usage{data(ExpressionControlData)}

\examples{
library(beadarray)
data(ExpressionControlData)
names(ExpressionControlData)
ExpressionControlData[["Humanv3"]][1:10,]

}

\keyword{datasets}

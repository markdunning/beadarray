\name{limmaDE}
\alias{limmaDE}

\title{
Differential expression using limma
}
\description{
Function to perform a standard limma analysis using a single command.
}
\usage{
limmaDE(summaryData, SampleGroup = NULL, DesignMatrix = NULL, makeWts = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{summaryData}{
ExpressionSetIllumina object
}
  \item{SampleGroup}{
vector to describe the sample groups to be compared
}
  \item{DesignMatrix}{
Optional design matrix
}
  \item{makeWts}{
if TRUE weights will be calculated for each array
}
  \item{\dots}{
Other arguments that lmFit can accept
}
}
\details{
All pairwise contrasts will be made
}
\value{
a limmaResults object
}
\references{

}
\author{
Mark Dunning
}

\examples{

if(require(beadarrayExampleData)){
  data(exampleSummaryData)
  rawdata <- channel(exampleSummaryData, "G")
  normdata <- normaliseIllumina(rawdata)
  limmaResults <- limmaDE(normdata, SampleGroup = pData(normdata)$SampleFac)
  limmaResults
}

}


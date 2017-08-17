\name{limmaDE}
\alias{limmaDE}

\title{
Differential expression using limma
}
\description{
Function to perform a standard limma analysis using a single command.
}
\usage{
limmaDE(summaryData, SampleGroup, DesignMatrix = NULL, makeWts = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{summaryData}{
ExpressionSetIllumina object
}
  \item{SampleGroup}{
Name of column in phenoData that will be used to construct sample groups for analysis
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
The function automates the steps used in a typical limma analysis. Firstly, the \link[limma]{lmFit} is used to fit a linear model from the specified column in \code{phenoData}. Array weights can be calculated (\link[limma]{arrayWeights}) and used in the fit. A contrast matrix of all possible contrasts is created (\link[limma]{makeContrasts}) and fitted (\link[limma]{makeContrasts}). The empirical Bayes shrinkage of variances is then applied (\link[limma]{eBayes}). The design matrix, contrast matrix and array weights can all be extracted for the resulting object.
}
\value{
a limmaResults object
}

\author{
Mark Dunning
}

\examples{

if(require(beadarrayExampleData)){
  data(exampleSummaryData)
  rawdata <- channel(exampleSummaryData, "G")
  normdata <- normaliseIllumina(rawdata)
  limmaResults <- limmaDE(normdata, SampleGroup = "SampleFac")
  limmaResults
  DesignMatrix(limmaResults)
}

}


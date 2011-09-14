\name{illuminaChannel-class}
\docType{class}
\alias{illuminaChannel-class}
\alias{greenChannel}

\title{Class "illuminaChannel" }

\description{
    A class to define how illumina bead-level data are summarized
}


\section{Slots/List Components}{
  Objects of this class contain the following slots
  \tabular{ll}{
  \code{transFun}:\tab function to transform the data from each array-section.\cr
  \code{outlierFun}:\tab A function for identifying outliers from a list of bead intensiites and associated ArrayAddressIDs .\cr
  \code{exprFun}:\tab A function for producing a single summary of expression level from a vector of bead-type intensities. e.g. mean .\cr
  \code{varFun}:\tab A function for producing a single summary of variability from a vector of bead-type intensities. e.g. sd \cr
  \code{name}: \tab Character vector that defines a name for the channel
  }
}

\details{

From beadarray version 2.0 onwards, users are allowed more flexibility in how to create summarized data from bead-level data. The \code{illuminaChannel} is a means of allowing this flexibility by definining how summarization will be performed on each array section in the bead-level data object. The three keys steps applied to each section are; 1) use a transform function to get the quantities to be summarized (one value per bead). The most common use-case would be to extract the Green channel intensities and possibly perform a log2 transformation. 2) remove any outliers from this list of values 3) split the values according to ArrayAddressIDs and apply the definied \code{exprFun} and \code{varFun} to the quantities belonging to each ArrayAddress.


}

\author{Mark Dunning}

\examples{

greenChannel

redChannel = new("illuminaChannel", redChannelTransform, illuminaOutlierMethod, mean, sd, "R")

logRatio = new("illuminaChannel", logRatioTransform, illuminaOutlierMethod, mean, sd, "M")

}




\seealso{\code{\link{summarize}}}
\keyword{classes}

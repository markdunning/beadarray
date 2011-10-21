\name{makeQCTable}
\alias{makeQCTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tabulate QC scores
}
\description{
Function to make a table of quality control scores for every section in a \code{beadLevelData} object. Either the annotation of the data needs to be specified, or a control profile data frame that lists ArrayAddress IDs and control types. The supplied summary functions are applied to each control type on each section.
}
\usage{
makeQCTable(BLData, transFun = logGreenChannelTransform, controlProfile = NULL, summaryFns = list(Mean = function(x) mean(x, na.rm=TRUE), Sd = function(x) sd(x, na.rm=TRUE)), channelSuffix = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{
a \code{beadLevelData} object
}
  \item{transFun}{
a function to be applied to the \code{beadLevelData} object prior to tabulation
}
  \item{controlProfile}{
an optional data frame that specifies ID and type of controls to be used to generate the table
}
  \item{summaryFns}{
list of functions to apply to each control type on every section
}
  \item{channelSuffix}{
optional character string to append to the column names of the resulting table
}
}
\details{

For each section in turn, the function groups together IDs of the same control type (e.g. housekeeping), and uses an \code{\link{lapply}} with the specified summary functions. A transformation function is applied to \code{BLData} prior to the summary, with the default being to take the log2 of the green channel.

If the annotation of the \code{beadLevelData} has been set by \code{\link{readIllumina}} or \code{\link{setAnnotation}} then the controlProfile data frame is calculated automatically and the \code{controlProfile} argument may be omitted.
}

\value{
A matrix with one row per section and one column for each combination of control type and summary function.
}

\examples{

if(require(beadarrayExampleData)){

data(exampleBLData)


qct = makeQCTable(exampleBLData)

qct

}

}

\author{
Mark Dunning
}



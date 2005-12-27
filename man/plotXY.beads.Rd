\name{plotXY.beads}
\alias{plotXY.beads}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{XY plots for two samples}
\description{
Function to display an XY plot of the average values of genes from two different samples or from two channels from the same array.

If two array indices are given as parameters then two arrays will be compared, otherwise the red and green channels from one array will be plotted
}
\usage{
plotXY.beads(BSData, array1, array2 = 0, log = FALSE, identify = FALSE, label = FALSE, fold = 2, highlight = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BSData}{BeadSummaryList object} 
  \item{array1}{numeric value for the first sample to be plotted}
  \item{array2}{numeric value for the second sample to be plotted}
  \item{log}{if TRUE all data will be log2-transformed}
  \item{identify}{if TRUE the user will be able to find additional information about each point by clicking on it}
  \item{label}{if TRUE the position of control types will be shown on the plot}
  \item{fold}{numeric value specifying what fold change lines to display on the plot}
  \item{highlight}{vector of numeric values that the user wishes to be highlighted}
}
\details{

A plot will be produced of the average values for each gene found in sample 1 against the averages found in sample 2. 

If identify is set to TRUE then clicking on each point will display information about the probeID, gene name and values from both samples of the point 
clicked

}
\value{
Plot will be produced on the current graphical device
}

\author{Mark Dunning}
\keyword{hplot}
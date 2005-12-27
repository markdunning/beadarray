\name{plotMA.beads}
\alias{plotMA.beads}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Show MA plots}
\description{
Function to show an MA plot between two specified samples given by different SAMs or the red and green channels of the same sample.

If two array indices are given then the MA plot will be between the two arrays, otherwise the two channels of one array will be plotted.

}
\usage{
plotMA.beads(BSData, array1, array2 = 0, identify = FALSE, label =
FALSE, highlight = NULL, log = TRUE, main = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BSData}{BeadSummaryList object}
  \item{array1}{numeric value to specify the first sample to be plotted}
  \item{array2}{numeric value to specify the second sample to be plotted (optional)}
  \item{identify}{if TRUE the user will be able to click on points on the plot and bring up useful information, such as probeID, gene information etc}
  \item{label}{if TRUE the control types for the experiment will be displayed on the plot}
  \item{highlight}{vector of numeric values specifying any probeIDs which the user wishes to be highlighted on the plot}
  \item{log}{if TRUE then a log2 transformation will be applied to the data before plotting}
  \item{main}{specifies a title to be written on the plot}
}
\details{



In the usual way, M is labelled on the y axis and represents the difference between two genes in different samples where as A is the average of the two samples.

Control information must be supplied in a controlTypes vector beforehand.
}
\value{
MA plot is displayed on current graphical device.
 
}

\author{Mark Dunning}

\keyword{hplot}
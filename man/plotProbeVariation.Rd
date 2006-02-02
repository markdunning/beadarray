\name{plotProbeVariation}
\alias{plotProbeVariation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot variation}
\description{
Plot the variation in the averages for a particular probe / bead type
}
\usage{
plotProbeVariation(BSData, probeID, log=TRUE, ylim=range(0,16), main=probeID,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BSData}{an BeadSummaryList containing bead summary data}
  \item{probeID}{character string specifying the probeID to be plotted}
  \item{log}{if TRUE the values will be transformed to be log2 scale before plotting}
  \item{ylim}{sets the range of the y axis on the plot}
  \item{main}{title for the plot}
  \item{\dots}{other graphical parameters to be passed}
}
\details{
The function plots the average (log2 or unlogged) value for the specified probeID against array index. This function is more suited to BeadChip data where the number of arrays
is relatively small (3 or 4). 

The plotOnSAM function can be used for the results of SAM experiments to relate average probe values to position on the SAM.

}

\value{
  Plot is produced on current graphical device
  
}

\author{Mark Dunning}

\keyword{methods}
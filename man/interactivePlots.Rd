\name{interactivePlots}
\alias{BeadChipSummary}
\alias{SAMSummary}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Interactive bead level plotting}
\description{
Plots using bead level data to discover spatial effects and artifacts on arrays

}
\usage{

SAMSummary(BLData, mode="outliers", missing_arrays=NULL, colour=TRUE, scale = NULL,...)
BeadChipSummary(BLData,  mode="outliers",colour=TRUE, scale = NULL,...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{BeadLevelList object} 
  \item{mode}{character string either "outliers", "fg", or "bg" specifying what to display on the plots}  
  \item{scale}{numeric value giving the amount by which to divide all numbers by (eg for log2 intensities this should be 16) to transform to range 0 - 1}
  \item{missing_arrays}{vector of numeric values specifying the index of any arrays that have been removed
    from the SAM.}
  \item{colour}{if TRUE the hexagons will be plotted in colour}
  \item{\dots}{other graphical parameters to plot that can be specified}

}
\details{

A plot will be displayed giving a summary of each array in the experiment on the left screen and initially a blank right hand side. The left hand side
is coloured according to the number of outliers found on the array or the mean intensity of the array (depending on the mode parameter). Clicking on tha particular 
array on the left will display a location plot of the outliers or an image plot on the right

See RNews article


}
\value{

Plot will be produced on the current graphical device

}

\author{Mark Dunning}
\keyword{hplot}
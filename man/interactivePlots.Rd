\name{interactivePlots}
\alias{BeadChipSummary}
\alias{SAMSummary}
\alias{screenSetup}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Interactive bead level plotting}
\description{
Plots using bead level data to discover spatial effects and artifacts on arrays

}
\usage{

SAMSummary(BLData, mode="outliers")
BeadChipSummary(BLData, mode="outliers")

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{BeadLevelList object} 
  \item{mode}{character string either "outliers", "fg", or "bg" specifying what to display on the plots}

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
\name{BLImagePlot}
\alias{BLImagePlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Image plot of BLData object}
\description{
Plots an image plot of a bead level data.  
}
\usage{
BLImagePlot(BLData, array = 1, nrow = 18, ncol = 2, low= NULL,
high = NULL, ncolors = 123)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{Object of class BLData}
  \item{array}{array to plot}
  \item{nrow}{How many rows should the array be divided into}
  \item{ncol}{How many columns should the array be divided into}
  \item{low}{Color for area of lowest intensity}
  \item{high}{color for area of highest intensity}
  \item{ncolors}{The number of color graduations between high and low}
}
\details{
Because of the high volume
of beadlevel data, this function works by mapping a grid (whos size is specified by the nrow and
ncol) arguments and averages the intensities of the beads within each
section of this grid.
}
\value{
Plot is produced on current graphical device.
}

\author{Mike Smith}

\keyword{hplot}
\name{imageplot}
\alias{imageplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{imageplot for BeadArrays object}
\description{
Plots an image plot of a bead level data object.  
}
\usage{
imageplot(BLData, array = 1, nrow = 18, ncol = 2, low= NULL,
high = NULL, ncolors = 123, whatToPlot="G",...)
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
  \item{whatToPlot}{string identifying which matrix in BLData to use in the imageplot}
  \item{\dots}{other graphical parameters to plot that can be specified}
}
\details{
Because of the high volume
of bead level data, this function works by mapping a grid (who's size is specified by the nrow and
ncol) arguments and averages the intensities of the beads within each
section of this grid.

The number of rows and columns may change the appearance of the plots. If the array is divided into too many sqaures it will difficult to detect changes. 
We recommend using nrow=200 and ncol=20 for whole genome arrays

An imageplot of the forgeround intensities is produced be default. Other matrices in BLData can be plotted by using the whatToPlot argument. eg change to
Gb to plot the background intensities.

}
\value{
Plot is produced on current graphical device.
}

\author{Mike Smith, Mark Dunning}

\keyword{hplot}
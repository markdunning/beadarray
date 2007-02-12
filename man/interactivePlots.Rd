\name{interactivePlots}
\alias{BeadChipSummary}
\alias{SAMSummary}

\title{Interactive bead level plotting}
\description{
Plots using bead level data to discover spatial effects and artifacts on arrays

}
\usage{

SAMSummary(BLData, mode = "outliers", whatToPlot = "G", log = TRUE,
               n = 3, missing_arrays = NULL, colour = TRUE,
               scale = NULL, low = "yellow", high="red", ...)
BeadChipSummary(BLData, mode = "outliers", whatToPlot = "G", log = TRUE,
               n = 3, colour = TRUE, scale = NULL, low = "yellow",
               high = "red", ...)
}

\arguments{
  \item{BLData}{BeadLevelList object}
  \item{mode}{character string either "outliers" or "intensities"
    specifying what to display on the plots}
  \item{whatToPlot}{character string specifying which intensities to plot.
    Possibilities are \code{"G"}, \code{"Gb"} for single channel data
    and \code{"G"}, \code{"Gb"}, \code{"R"} and \code{"Rb"} for
    two-colour data}
  \item{log}{if TRUE log2 intensities of each bead are used to
    find outliers}
  \item{n}{Specify a cut-off for outliers as n median
    absolute deviations from the median. The default value is 3}
  \item{missing_arrays}{vector of numeric values specifying the index
    of any arrays that have been removed from the SAM.}
  \item{colour}{if TRUE the hexagons will be plotted in colour}
  \item{scale}{numeric value giving the amount by which to divide all
    numbers by (eg for log2 intensities this should be 16) to transform
    to range 0 - 1}
   \item{low}{Color for area of lowest intensity}
  \item{high}{color for area of highest intensity}
  \item{\dots}{other graphical parameters to plot that can be specified}
}
\details{
  A plot will be displayed giving a summary of each array in the
  experiment on the left screen and initially a blank right hand side.
  The left hand side is coloured according to the number of outliers
  found on the array or the mean intensity of the array (depending on
  the mode parameter). Clicking on tha particular array on the left
  will display a location plot of the outliers or an image plot on the
  right See RNews article
}
\value{
  Plot will be produced on the current graphical device
}

\author{Mark Dunning}
\keyword{hplot}
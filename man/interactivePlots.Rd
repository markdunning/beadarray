\name{interactivePlots}
\alias{BeadChipSummary}
\alias{SAMSummary}

\title{Interactive bead-level plotting}
\description{
Generates spatial plots using bead-level data to discover artefacts on strips/arrays.
}
\usage{

SAMSummary(BLData, mode = "outliers", whatToPlot = "G",  samID = NULL, 
               log = TRUE, n = 3, colour = TRUE,
               scale = NULL, low = "yellow", high="red", ...)
BeadChipSummary(BLData, mode = "outliers", whatToPlot = "G", chipID = NULL, 
               stripsPerChip = 12, log = TRUE, n = 3, colour = TRUE, 
               scale = NULL, low = "yellow", high = "red", ...)
}

\arguments{
  \item{BLData}{\code{BeadLevelList} object}
  \item{mode}{character string either "outliers" or "intensities"
    specifying what to display on the plots}
  \item{whatToPlot}{character string specifying which intensities to plot.
    Possibilities are \code{"G"}, \code{"Gb"} for single channel data
    and \code{"G"}, \code{"Gb"}, \code{"R"} and \code{"Rb"} for
    two-colour data}
  \item{samID}{character string specifying which SAM to plot.  If \code{NULL}, data from the first SAM is plotted.}
  \item{chipID}{character string specifying which BeadChip to plot.  If \code{NULL}, data from the first BeadChip is plotted.}
 \item{stripsPerChip}{integer specifying number of strips on BeadChip (8 or 12)}
  \item{log}{if TRUE log2 intensities of each bead are used to
    find outliers}
  \item{n}{numeric value specifying the number of median absolute deviations (MADs)
    from the median to use as a cut-off for outliers.  The default value is 3}
  \item{colour}{if TRUE the hexagons will be plotted in colour}
  \item{scale}{numeric value giving the amount by which to divide all
    numbers by (eg for log2 intensities this should be 16) to transform
    to range 0 - 1}
   \item{low}{colour to use for lowest intensity}
  \item{high}{colour to use for highest intensity}
  \item{\dots}{other parameters to \code{imageplot} that can be specified}
}
\details{
  A plot will be displayed giving a summary of each array in the
  experiment on the left screen and initially a blank right hand side.
  The left hand side is coloured according to the number of outliers
  found on the array or the mean intensity of the array (depending on
  the \code{mode} parameter). Clicking on tha particular array on the left
  will display a location plot of the outliers or an image plot on the
  right See RNews article
}
\value{
  A plot is produced on the current graphical device
}

\author{Mark Dunning}
\keyword{hplot}
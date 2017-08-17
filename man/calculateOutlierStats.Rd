\name{calculateOutlierStats}
\alias{calculateOutlierStats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Outlier distribution stats
}
\description{
Function that determines the outlier beads on an array and how they are distributed among the segments

}
\usage{
calculateOutlierStats(BLData, array = array, transFun = logGreenChannelTransform, outlierFun = illuminaOutlierMethod, n=3, useLocs = TRUE, nSegments = 9,...)
}

\arguments{
  \item{BLData}{
a \code{\link{beadLevelData-class}} object
}
  \item{array}{
the number of the array of interest
}
  \item{transFun}{
how the section data is to be transformed prior to calculating outliers
}
  \item{outlierFun}{
a function for calculating outliers
}
  \item{n}{
an indicator of how extreme an observation must be (e.g. how many MADs from the median), to be passed to the function that will identify outliers}
  \item{useLocs}{
use locs and sdf information (if available) to determine section layout
}
  \item{nSegments}{
manually set how many segments the section is divided into
}
  \item{...}{Additional arguments to be passed to \code{outmeth}.}
}
\details{

A section of an expression BeadChip (e.g. the Humanv3 or HumanHT-12) is made up of 9 physcially-separate segments. A useful QA check is to see how the outliers are distributed among these segments. Outliers are beads that have outlying intensities according to some rule that the user can specify. The default (as used by Illumina) is to exclude beads that are more than 3 median absolute deviations from the median. Once outliers are determined, the coordinates for these outliers are binned into segments by assuming that the segments are evenly spaced across the section surface.

Note that sections from Sentrix Array Matrix do not have segments, so the results may not be informative

}


\value{
vector with the percentage of beads found in each segment that were determined to be outliers

}


\author{
Mark Dunning
}

\examples{

if(require(beadarrayExampleData)){

data(exampleBLData)

##Artificial example, there are no segments on this type of BeadArray

calculateOutlierStats(exampleBLData, array=1, nSegments=10, useLocs=FALSE)

calculateOutlierStats(exampleBLData, array=2, nSegments=10, useLocs=FALSE)


}


}


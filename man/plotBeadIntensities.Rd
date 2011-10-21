\name{beadIntensityPlots}
\alias{plotBeadIntensities}
\alias{genericBeadIntensityPlot}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting the intensities of selected beads on a section
}
\description{
The function will plot the intensities of selected beads on a specified array
}
\usage{
plotBeadIntensities(BLData, array = 1, BeadIDs, transFun = logGreenChannelTransform, cols = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{
a \code{beadLevelData} object
}
  \item{array}{
numeric specifying which array to plot the intensities from
}
  \item{BeadIDs}{
what ArrayAddress IDs to be plotted
}
  \item{transFun}{
function specifying what transformation to be applied to the \code{beadLevelData} prior to plotting 
}
  \item{cols}{
a vector of colours to be used to plot each ID. If NULL the \code{rainbow} function is used to generate colours.

}
  \item{\dots}{
other argument that may be passed along to plot.
}
}
\details{
The function will take all data from the specified section, apply the transformation (the default is to do log2) and then find the subset of beads that have the specified ID. These IDs should match the numeric ArrayAddress IDs that are stored in the beadLevelData object. 

}
\value{
Plot is produced on current graphical device.

}
\author{
Mark Dunning
}

\examples{

if(require(beadarrayExampleData)){


data(exampleBLData)

randIDs = sample(getBeadData(exampleBLData, array=1, what="ProbeID"),10)

plotBeadIntensities(exampleBLData, array=1, BeadIDs = randIDs)

}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplots}


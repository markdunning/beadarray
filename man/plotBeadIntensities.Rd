\name{plotBeadIntensities}
\alias{plotBeadIntensities}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot bead intensities}
\description{
Function to plot the intensities of all beads of a particular type on an array against their distance from the centre of the array. Any outliers 
for the bead type will be highlighted.

}
\usage{
plotBeadIntensities(BLData, probe, array, n=3, log=FALSE, label=FALSE, ylim
= NULL, main = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{BeadLevelList object}
  \item{probe}{numeric value giving the ProbeID for the bead type we are interested in}
  \item{array}{numeric value giving the array we are interested in}
  \item{n}{outlier values for the bead type are beads with intensity more than 'n' times away from the mean}
  \item{log}{if TRUE log2 intensities of each bead are used to find outliers}
  \item{label}{if TRUE all beads will be given a numeric ID on the plot}
  \item{ylim}{vector providing the upper and lower bounds for the
    y-axis in the plotting window}
\item{main}{specify a title for the plot}
\item{\dots}{further graphical parameters to \code{plot}}
  
}
\details{
We find the intensities for each bead of the specified bead type by simply finding the value of BLData\$R for each bead in the array whose BLData\$ProbeID 
matches the specified ProbeID. The distance from the centre of the hexagon can be found in a similar fashion by looking at BLData\$x and BLData\$y. 

The outliers for the bead type are beads which have intensity more than 'n' times away from the standard deviation using logged or un-logged bead intensities 
(specified using 'log' argument).

The mean for the bead type and lines 3 standard deviations from the mean are shown on the plot.

}
\value{
Plot is produced on current graphical device
}

\author{Mark Dunning}

\examples{

data(BLData)

#By default, outliers found by Illumina method are highlighted

plotBeadIntensities(BLData, 2, 1)

plotBeadIntensities(BLData, 2, 1, log=TRUE)

plotBeadIntensities(BLData, 2, 1, log=TRUE, n=5)

#Increasing n selects only the most extreme outliers

plotBeadIntensities(BLData, 4279, 1, log=TRUE, n=3)

plotBeadIntensities(BLData, 4279, 1, log=TRUE, n=4.5)

plotBeadIntensities(BLData, 4279, 1, log=TRUE, n=6)



}
\keyword{hplot}
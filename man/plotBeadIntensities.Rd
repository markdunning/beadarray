\name{plotBeadIntensities}
\alias{plotBeadIntensities}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot bead intensities}
\description{
Function to plot the intensities of all beads of a particular type on an array against their distance from the centre of the array. Any outliers 
for the bead type will be highlighted.

}
\usage{
plotBeadIntensities(BLData, ProbeIDs, arrays, log=FALSE, n=3,ProbeCols=NULL,ylim=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{BeadLevelList object}
  \item{ProbeIDs}{numeric value giving the ProbeID for the bead type we are interested in}
  \item{arrays}{numeric value giving the array we are interested in}
  \item{log}{if TRUE log2 intensities of each bead are used to find outliers}
  \item{n}{outlier values for the bead type are beads with intensity more than 'n' times away from the mean}
  \item{ProbeCols}{vector of colour names to assign to each bead type}
  \item{ylim}{vector providing the upper and lower bounds for the
    y-axis in the plotting window}
\item{\dots}{further graphical parameters to \code{plot}}
  
}
\details{
We find the intensities for each bead of the specified bead type by simply finding the value of BLData\$G for each bead in the array whose BLData\$ProbeID 
matches the specified ProbeID. The distance from the centre of the hexagon can be found in a similar fashion by looking at BLData\$GrnX and BLData\$GrnY. 

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

ids = intersect(BLData@ProbeID[,1], BLData@ProbeID[,3])[1:10]

plotBeadIntensities(BLData, ProbeIDs=ids[1], arrays=1, ylim=range(9,11))

plotBeadIntensities(BLData, ProbeIDs=ids[1], arrays=c(1,3), ylim=range(9,11))

plotBeadIntensities(BLData, ProbeIDs=ids, arrays=c(1,3), ylim=range(9,11))


}
\keyword{hplot}
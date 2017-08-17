\name{plotBeadLocations}
\alias{plotBeadLocations}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot bead locations
}
\description{
Can plot where specified beads, or bead types were located on the array surface
}
\usage{
plotBeadLocations(BLData, ProbeIDs = NULL, BeadIDs = NULL, array = 1, SAM = FALSE, xCol = "GrnX", yCol = "GrnY", xlab = "x-coordinate", ylab = "y-coordinate", horizontal = TRUE, main = paste("Bead", ProbeIDs, "locations"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{
a \code{beadLevelData} object
}
  \item{ProbeIDs}{
a list of ArrayAddress IDs to plot
}
  \item{BeadIDs}{
a list of beads (rows in the beadLevelData object) to plot
}
  \item{array}{
the number of the section to plot
}
  \item{SAM}{
if TRUE the array is treated as a Sentrix Array Matrix (hexagonal)
}
  \item{xCol}{
column name for the x coordinates
}
  \item{yCol}{
column name for the x coordinates
}
  \item{xlab}{
optional label for the x axis
}
  \item{ylab}{
optional label for the y axis
}
  \item{horizontal}{
if TRUE the longest edge of the array surface will be plotted on the x axis
}
  \item{main}{
an optional title for the plot
}
  \item{\dots}{
any arguments to be passed to plots
}
}
\value{
plot to current graphical device
}
\author{
Mark Dunning
}

\examples{

if(require(beadarrayExampleData)){

data(exampleBLData)

##Plot location of first 100 beads as they are listed in beadLevelData object

plotBeadLocations(exampleBLData, BeadIDs = 1:100, array=1, horizontal = FALSE)



}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

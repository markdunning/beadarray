\name{poscontPlot}
\alias{poscontPlot}
\title{
Plot the positive controls 
}
\description{
Function for retrieving and plotting the biotin and housekeeping controls for an expression array. We know these controls should show high signal and are therefore useful for QA purposes. The housekeeping control targets a bead-type believed to be universally expressed whereas the biotin control targets the biotin used for staining. 
}
\usage{
poscontPlot(BLData, array = 1, transFun = logGreenChannelTransform, positiveControlTags = c("housekeeping", "biotin"), colList = c("red", "blue"), controlProfile = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{
a \code{beadLevelData} object
}
  \item{array}{
The section to be plotted
}
  \item{transFun}{
What transformation function to be applied prior to plotting
}
  \item{positiveControlTags}{
What identifiers to be used as positive controls
}
  \item{colList}{
vector of colours to be used to each positive control
}
  \item{controlProfile}{
an optional data frame with columns defining the ArrayAddress IDs and control\-type for all controls on the platform.
}
  \item{\dots}{
other arguments to plot
}
}
\details{

Function for plotting the observed intensites for all replicates of the specified control probes on a given array\-section. The identity of the control probes can be specified by passing a ControlProfile data frame, with the first column being a vector of ArrayAddress IDs and the second column being a corresponding set of characters tags. The beads to be plotted are found by matching the \code{positiveControlTags} argument to these character tags. 
Users with expression data can have the ControlProfile data frame defining automatically within the function, provided the annotation of the \code{beadLevelData} object has been defined by \code{readIllumina} or \code{setAnnotation}.

}
\value{
Plot to current graphical device

}
\references{\url{www.illumina.com/downloads/GX\_QualityControl\_TechNote.pdf}}

\author{
Mark Dunning
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

###Load the example beadLevelData and associated controlProfile

if(require(beadarrayExampleData)){


data(exampleBLData)

poscontPlot(exampleBLData, array=1)

poscontPlot(exampleBLData, array=2)

}

}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplots}


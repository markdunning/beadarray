\name{maplots}
\alias{maplots}

\title{
Plot for pairwise comparisons of array intensities
}
\description{
Function to construct the classic MA plots from a dataset
}
\usage{
maplots(data, sampleFactor = NULL, do.log = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
An ExpressionSetIllumina object
}
  \item{sampleFactor}{
If left NULL, all arrays will be compared to a reference defined by the mean intensities across all arrays. Otherwise, all pairwise comparisons will be made within each sample group separately.
}
  \item{do.log}{
if TRUE, a log2 transform will be performed prior to plotting.
}
}
\details{
Densities are computed using the 'hexbin' package
}
\value{
A list of ggplot2 objects is returned.

}

\author{
Mark Dunning
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

if(require(beadarrayExampleData)){
data(exampleSummaryData)
rawdata <- channel(exampleSummaryData, "G")
 mas <- maplots(rawdata,do.log=FALSE)
mas[[1]]
##Smoothed lines can be added
mas[[1]] + geom_smooth(col="red")
##Added lines on the y axis
mas[[1]] + geom_hline(yintercept=c(-1.5,1.5),col="red")
##Changing the color scale
mas[[1]] + scale_fill_gradient2(low="yellow",mid="orange",high="red")

}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

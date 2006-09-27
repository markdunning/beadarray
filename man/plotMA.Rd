\name{plotMA}
\alias{plotMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Show MA plots}
\description{
Function to show an MA plot between two specified samples given by different SAMs or the red and green channels of the same sample.

If two array indices are given then the MA plot will be between the two arrays, otherwise the two channels of one array will be plotted.

}
\usage{
plotMA(exprs, array1, array2=0, genesToLabel=NULL, labelCol="red", foldLine=2,log=TRUE, labelpch=16,identify=FALSE, ma.ylim=2,sampleSize=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exprs}{expression matrix from an ExpressionSet or eSet}
  \item{array1}{numeric value to specify the first sample to be plotted}
  \item{array2}{numeric value to specify the second sample to be plotted (optional)}
  \item{genesToLabel}{vector of gene names to highlight on the plot. These must match the rownames of exprs}
  \item{labelCol}{vector of characters to plot for the selected genes}
  \item{log}{if TRUE the data will be log-transformed before plotting}
  \item{labelpch}{if TRUE then a log2 transformation will be applied to the data before plotting}
  \item{foldLine}{a numeric value defining where to draw horizontal fold change lines on the plot}
  \item{identify}{if TRUE, clicking on points on the plot will bring up extra information}
  \item{sampleSize}{numeric specifiying the size of a random sample to be taken rather than plotting every gene}
  \item{ma.ylim}{numeric value specifying the range of the plot (from -ma.ylim to ma.ylim)}
  \item{\dots}{other graphical parameters to plot that can be specified}
}
\details{


In the usual way, M is labelled on the y axis and represents the difference between two genes in different samples where as A is the average of the two samples.

Control information must be supplied in a controlTypes vector beforehand.
}
\value{
MA plot is displayed on current graphical device.
 
}

\author{Mark Dunning}

\keyword{hplot}
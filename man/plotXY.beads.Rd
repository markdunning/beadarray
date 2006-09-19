\name{plotXY}
\alias{plotXY}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{XY plots for two samples}
\description{
Function to display an XY plot of the average values of genes from two different samples or from two channels from the same array.

If two array indices are given as parameters then two arrays will be compared, otherwise the red and green channels from one array will be plotted
}
\usage{
plotXY(exprs,array1, array2=0, genesToLabel=NULL, labelCol="red", labelpch=16,identify=FALSE, foldLine=2,sampleSize=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exprs}{expression matrix from an ExpressionSet or eSet object} 
  \item{array1}{numeric value for the first sample to be plotted}
  \item{array2}{numeric value for the second sample to be plotted}
  \item{genesToLabel}{vector of gene names to highlight on the plot. These must match the rownames of exprs}
  \item{labelCol}{vector of characters to plot for the selected genes}
  \item{labelpch}{if TRUE then a log2 transformation will be applied to the data before plotting}
  \item{foldLine}{a numeric value defining where to draw horizontal fold change lines on the plot}
  \item{identify}{numeric value specifying what fold change lines to display on the plot}
  \item{sampleSize}{numeric specifiying the size of a random sample to be taken rather than plotting every gene}
  \item{\dots}{other graphical parameters that can be sent to plot}
}
\details{

A plot will be produced of the average values for each gene found in sample 1 against the averages found in sample 2. 

If identify is set to TRUE then clicking on each point will display information about the ProbeID, gene name and values from both samples of the point 
clicked

}
\value{
Plot will be produced on the current graphical device
}

\author{Mark Dunning}
\keyword{hplot}
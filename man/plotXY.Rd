\name{plotXY}
\alias{plotXY}
\title{XY plots for two samples}
\description{
Function to display an XY plot of the average values of genes from two different samples or from two channels from the same array.

If two array indices are given as parameters then two arrays will be compared, otherwise the red and green channels from one array will be plotted
}
\usage{
plotXY(exprs,array1=1, array2=2, genesToLabel=NULL, labelCol="red", log=TRUE,labelpch=16, foldLine=2,sampleSize=NULL,...)
}
\arguments{
  \item{exprs}{expression matrix from an ExpressionSet or eSet object} 
  \item{array1}{numeric value for the first sample to be plotted}
  \item{array2}{numeric value for the second sample to be plotted}
  \item{genesToLabel}{vector of genes to highlight on the plot. These must match the rownames of exprs}
  \item{labelCol}{plotting colours for highlighted genes}
  \item{log}{if TRUE the data will be log-transformed before plotting}
  \item{labelpch}{plotting characters for highlighted genes}
  \item{foldLine}{a numeric value defining where to fold change lines on the plot}
  \item{sampleSize}{The number of genes to plot.  Default is NULL, which plots every gene.}
  \item{\dots}{other graphical parameters to plot that can be specified}
}
\details{

Plots the \code{array1} intensity versus the \code{array2} intensity for each gene. 

}
\value{
Plot will be produced on the current graphical device
}
\author{Mark Dunning}

\examples{
data(BSData)

plotXY(exprs(BSData), array1=1, array2=2)
}
\keyword{hplot}

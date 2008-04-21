\name{plotXY}
\alias{plotXY}
\title{XY plots for two samples}
\description{
Function which produces an XY plot of the intensities from two specified arrays.
}
\usage{
plotXY(exprs,array1=1, array2=2, genesToLabel=NULL, labelCol="red", log=TRUE,labelpch=16, foldLine=2,sampleSize=NULL,...)
}
\arguments{
  \item{exprs}{a matrix of expression values} 
  \item{array1}{integer specifying the first array to plot}
  \item{array2}{integer specifying the second array to plot}
  \item{genesToLabel}{vector of genes to highlight on the plot. 
    These must match the rownames of \code{exprs}.}
  \item{labelCol}{plotting colours for highlighted genes}
  \item{log}{if TRUE the data will be log-transformed before plotting}
  \item{labelpch}{plotting characters for highlighted genes}
  \item{foldLine}{a numeric value defining where to fold change lines on the plot}
  \item{sampleSize}{The number of genes to plot.  Default is NULL, which plots every gene.}
  \item{\dots}{other graphical parameters to plot that can be specified}
}
\details{
Plots the \code{array1} intensities versus the \code{array2} intensities for all probes. 
}
\value{
A smoothed scatter plot is displayed on the current graphical device
}
\author{Mark Dunning}

\examples{
data(BSData)

plotXY(exprs(BSData), array1=1, array2=2)
}
\keyword{hplot}

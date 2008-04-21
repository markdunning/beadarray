\name{plotMA}
\alias{plotMA}
\title{Show MA plots}
\description{
Function which produces an MA plot between two specified arrays.
}
\usage{
plotMA(exprs, array1=1, array2=2, genesToLabel=NULL, labelCol="red", foldLine=2, log=TRUE, labelpch=16, ma.ylim=2, sampleSize=NULL,...)
}
\arguments{
  \item{exprs}{a matrix of expression values}
  \item{array1}{integer specifying the first array to plot}
  \item{array2}{integer specifying the second array to plot}
  \item{genesToLabel}{vector of genes to highlight on the plot. 
     These must match the rownames of \code{exprs}.}
  \item{labelCol}{plotting colours for highlighted genes}
  \item{foldLine}{a numeric value defining where to draw horizontal fold-change lines on the plot}
  \item{log}{if TRUE the data will be log-transformed before plotting}
  \item{labelpch}{plotting characters for highlighted genes}
  \item{ma.ylim}{numeric value specifying the range of the plot (from -ma.ylim to ma.ylim)}
  \item{sampleSize}{The number of genes to plot.  Default is NULL, which plots every gene.}
  \item{\dots}{other graphical parameters to plot that can be specified}
}
\details{

The log2 difference in intensity (M-value, log-ratio) are plotted against the log2 average intensity (A-value) for each probe for the two arrays selected.
}
\value{
A smoothed MA scatter plot is displayed on the current graphical device.
}

\author{Mark Dunning}

\examples{
data(BSData)

plotMA(exprs(BSData), array1=1, array2=2)
}
\keyword{hplot}

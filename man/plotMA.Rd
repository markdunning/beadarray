\name{plotMA}
\alias{plotMA}
\title{Show MA plots}
\description{
Function to show an MA plot between two specified samples from  different BeadArrays.
}
\usage{
plotMA(exprs, array1=1, array2=2, genesToLabel=NULL, labelCol="red", foldLine=2, log=TRUE, labelpch=16, ma.ylim=2, sampleSize=NULL,...)
}
\arguments{
  \item{exprs}{expression matrix from an ExpressionSet or eSet}
  \item{array1}{numeric value to specify the first sample to be plotted}
  \item{array2}{numeric value to specify the second sample to be plotted (optional)}
  \item{genesToLabel}{vector of genes to highlight on the plot. These must match the rownames of exprs}
  \item{labelCol}{plotting colours for highlighted genes}
  \item{foldLine}{a numeric value defining where to draw horizontal fold change lines on the plot}
  \item{log}{if TRUE the data will be log-transformed before plotting}
  \item{labelpch}{plotting characters for highlighted genes}
  \item{ma.ylim}{numeric value specifying the range of the plot (from -ma.ylim to ma.ylim)}
  \item{sampleSize}{The number of genes to plot.  Default is NULL, which plots every gene.}
  \item{\dots}{other graphical parameters to plot that can be specified}
}
\details{

On the y axis, the log-2 difference in intensity (M-values) between the two samples is shown and along the x axis, the log-2 average intensity (A-values) are plotted.
}
\value{
MA plot is displayed on current graphical device.
 
}

\author{Mark Dunning}

\examples{
data(BSData)

plotMA(exprs(BSData), array1=1, array2=2)
}
\keyword{hplot}

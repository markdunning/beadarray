\name{plotMAXY}
\alias{plotMAXY}
\alias{plotMA}
\alias{plotXY}
\title{Scatter plots and MA-plots for all specified arrays}
\description{

Produces smoothed scatter plots of M versus A and X versus Y for all pairwise comparisons from a set of arrays.
}
\usage{
plotMAXY(exprs, arrays, log = TRUE, genesToLabel=NULL,
         labels=colnames(exprs)[arrays],labelCol="red",
         labelpch=16,foldLine=2,sampleSize=NULL,...)
}
\arguments{
  \item{exprs}{a matrix of expression values}
  \item{arrays}{integer vector giving the indices of the arrays 
    (columns of \code{exprs}) to plot}
  \item{log}{if TRUE then all values will be log2-transformed before plotting}
  \item{genesToLabel}{vector of genes to highlight on the plot.  
    These must match the rownames of \code{exprs}.}
  \item{labels}{vector of array names to display on the plot}
  \item{labelCol}{plotting colours for highlighted genes}
  \item{labelpch}{plotting characters for highlighted genes}
  \item{foldLine}{a numeric value defining where to draw horizontal fold change lines on the plot}
  \item{sampleSize}{The number of genes to plot.  Default is NULL, which plots every gene}
  \item{\dots}{other graphical parameters to be passed}
}

\details{
This graphical tool shows differences that exist between two arrays and can be used to highlight biases between arrays as well as highlighting genes which are differentially expressed.  For each bead type, we calculate the average (log2) intensity and difference in intensity (log2-ratio) for each pair of arrays.

In the lower-left section of the plot we see XY plots of the intensities for all pairwise comparisons between the arrays and in the upper right we have pairwise MA plots. Going down the first column we observe XY plots of array 1 against array 2 and array 1 against array 3 etc. Similarly, in the upper-right corner we can observe pairwise MA plots.
}
\author{Mark Dunning}

\examples{

if(require(beadarrayExampleData)){

data(exampleSummaryData)

plotMAXY(exprs(channel(exampleSummaryData, "G")), arrays=1:3, log=FALSE)

}

}
\keyword{hplot}

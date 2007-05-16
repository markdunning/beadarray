\name{plotMAXY}
\alias{plotMAXY}

\title{Scatter plots and MA-plots for all specified arrays}
\description{

Produces scatter and MA-plots for all pairwise comparisons that can be made from a set of arrays.
}
\usage{
plotMAXY(exprs, arrays, log = TRUE, genesToLabel=NULL,labels=NULL,labelCol="red", labelpch=16,foldLine=2,sampleSize=NULL,...)
}
\arguments{
  \item{exprs}{an expression matrix from an ExpressionSet or eSet}
  \item{arrays}{numeric vector giving the indices of the arrays in BSData we want to plot}  
  \item{log}{if TRUE then all values will be log2-transformed before plotting}
  \item{genesToLabel}{vector of genes to highlight on the plot. These must match the rownames of exprs}
  \item{labels}{vector of array names to display on the plot}
  \item{labelCol}{plotting colours for highlighted genes}
  \item{labelpch}{plotting characters for highlighted genes}
  \item{foldLine}{a numeric value defining where to draw horizontal fold change lines on the plot}
  \item{sampleSize}{The number of genes to plot.  Default is NULL, which plots every gene}
  \item{\dots}{other graphical parameters to be passed}
}

\details{
This graphical tool shows differences that exist between two arrays and can be used to highlight biases between arrays as well as highlighting genes which are differentially expressed.  For each bead type, we take the averaged intensity calculated from two different arrays. On an XY plot these values are plotted against each other and on the MA plot, log2-ratios is plotted against the average intensities.

In the lower-left section of the plot we see XY plots for all pairwise comparisons between the arrays and in the upper right we have pairwise MA plots. Going down the first column we observe XY plots of array 1 against array 2 and array 1 against array 3 etc. Similarly, in the upper-right corner we can observe pairwise MA plots.
}
\author{Mark Dunning}

\examples{
data(BSData)

vec = c(1,2,3)
labs = colnames(exprs(BSData))
plotMAXY(exprs(BSData), vec, labels=labs)

}
\keyword{hplot}

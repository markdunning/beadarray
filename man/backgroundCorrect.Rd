\name{backgroundCorrect}
\alias{backgroundCorrect}
\alias{backgroundCorrect,BeadLevelList-method}
\alias{backgroundCorrect,RGList-method}
\title{Background correct a BeadLevelList object}
\description{
  Adapted from the 'limma' backgroundCorrect function to correct the
  foreground intensities of a \code{BeadLevelList} object using the background values.
}
\usage{
backgroundCorrect(object, method = "subtract", offset = 0, verbose = FALSE)
}
\arguments{
  \item{object}{a BeadLevelList object}
  \item{method}{character string specifying correction method. Possible
          values are '"none"', '"subtract"', '"half"', '"minimum"',
	  '"edwards"', '"normexp"', '"rma"')}
  \item{offset}{numeric value to add to the intensities}
  \item{verbose}{logical. Used when method = '"normexp"'.  If TRUE,
    the parameters estimated by the model are output.}	
}
\details{
  Below is an excerpt from the 'limma' backgroundCorrect man page:

  If 'method="none"' then the corrected intensities are equal to the
  foreground intensities, i.e., the background intensities are
  treated as zero. If 'method="subtract"' then this function simply
  subtracts the background intensities from the foreground
  intensities which is the usual background correction method.

  The remaining methods are all designed to produce positive
  corrected intensities. If 'method="half"' then any intensity which
  is less than 0.5 after background subtraction is reset to be equal
  to 0.5. If 'method="minimum"' then any intensity which is zero or
  negative after background subtraction is set equal to half the
  minimum of the positive corrected intensities for that array.
  If 'method="edwards"' the method of Edwards (2003) is used.
  If 'method="normexp"' or '"rma"', a normal-exponential convolution
  model is fitted to the intensities, using different estimation
  procedures.  See Smyth (2005) for further details on normexp.

  The 'offset' can be used to add a constant to the intensities
  before log-transforming, so that the log-ratios are shrunk towards
  zero at the lower intensities. This may eliminate or reverse the
  usual 'fanning' of log-ratios at low intensities associated with
  local background subtraction.
  
  End of excerpt.
  
  As a result of both having identical function names this function can
  conflict with the
  \code{\link[limma:backgroundcorrect]{backgroundCorrect}} method in
  'limma'.  If both packages are loaded, the function from whichever
  package was loaded last takes precedence.  If the 'beadarray'
  backgroundCorrect() function is masking that from 'limma', one can
  directly call the 'limma' method using the command
  "limma::backgroundCorrect()".  Alternatively, one can detach the
  'beadarray' package using "detach(package:beadarray)".  Similar
  techniques can be used if 'limma' is masking the 'beadarray' method. 
}
\value{
  A \code{\link{BeadLevelList}} object in which the 'G' (and 'R', if present) intensities for each array are background corrected.  Note that the 'Gb' (and 'Rb' intensities) are not removed. 
}
\author{Mark Dunning and Mike Smith based on the limma function}
\references{
  Edwards, D. E. (2003). Non-linear normalization and background
  correction in one-channel cDNA microarray studies, Bioinformatics,
  19, 825-833.

  Smyth, G. K. (2005). Limma: linear models for microarray data. In:
  Bioinformatics and Computational Biology Solutions using R and
  Bioconductor, R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W.
  Huber (eds.), Springer, New York, pages 397-420.
}    
\examples{
data(BLData)

#default is to simply subtract Rb from R
BLData.bc = backgroundCorrect(BLData)

#Use 'minimum' method to stop negative values appearing
BLData.min = backgroundCorrect(BLData, method="minimum")
}
\keyword{methods}

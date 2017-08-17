\name{readIdatFiles}
\alias{readIdatFiles}
\title{Read BeadScan gene expression output}
\description{
Function to read IDAT files into beadarray.
}
\usage{
readIdatFiles(idatFiles = NULL)

}
\arguments{
  \item{idatFiles}{A vector of file paths to the idat files that are to be read.}
  }

\details{
  This function allows IDAT files, produced during the BeadArray scanning, 
  to be read directly into \code{beadarray}.  This removes the requirement to
  use Illumina's GenomeStudio software to convert such files into an ASCII format
  before they can be processed.
  
  The expection is that IDAT files from expression BeadArrays will be read using
  this function.  However, this is not currently checked for explicitly and it is
  possible that files from other platforms may be read successfully.

  If files from a two-colour array are to be read, they should all be listed 
  in the \code{idatFiles} argument.  The function will determine two-colour data
  is present and act accordingly.  If unequal numbers of green and red files are
  supplied, or the names of the array sections do not match, the function will 
  fail.
  
  IDAT files often contain probes that are used as internal controls and are not
  annotated by Illumina.  Such probes are not included in the output of this
  function.
}

\value{
  An \code{ExpressionSetIllumina} object.
}

\author{Mike Smith}

\keyword{IO}

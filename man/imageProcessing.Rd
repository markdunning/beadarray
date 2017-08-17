\name{imageProcessing}
\alias{illuminaForeground}
\alias{illuminaBackground}
\alias{medianBackground}
\alias{illuminaSharpen}

\title{Image processing functions}

\description{
Functions for obtaining bead intensity values from raw tiff images.  The three commands with the illumina prefix attempt to emulate the image processing implemented by Illumina.  The medianBackground function implements a more robust background calculation recommended by Smith et al.
}

\usage{
illuminaForeground(pixelMatrix, beadCoords)
illuminaBackground(pixelMatrix, beadCoords)
illuminaSharpen(pixelMatrix)

medianBackground(pixelMatrix, beadCoords)
}

\arguments{
  \item{pixelMatrix}{A matrix storing the individual pixel values of an image.  Intended to be created by \code{\link{readTIFF}}, although any matrix can be passed as input.}
  \item{beadCoords}{Two column matrix with each row containing a pair of coordinates representing a bead centre.}
}

\value{
illuminaForeground, illuminaBackground and medianBackground return a vector of intensity values, with one entry for every row in the beadCoords argument.  Any pairs of coordinates that fall outside the dimensions of the image return NA.

illuminaSharpen returns a matrix with the same dimensions as the pixelMatrix argument.
}

\author{Mike Smith}

\references{
Smith ML, Dunning MJ, Tavare S, Lynch AG. Identification and correction of previously unreported spatial phenomena using raw Illumina BeadArray data. \emph{BMC Bioinformatics} (2010) 11:208   
}

\keyword{manip}

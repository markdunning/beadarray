\name{readTIFF}
\alias{readTIFF}
\title{Read the Illumina tiff images}
\description{
  Reads Illumina tiff images and produces a matrix of pixel values.
}
\usage{
readTIFF(fileName, verbose = FALSE)
}

\arguments{
  \item{fileName}{String specifying the name of the tiff image to be read.}
  \item{verbose}{If TRUE then details from the header of the tiff are
  printed as it is read.  These include things like the byte order, the
  number of pixels in the image, the number of tags in the header etc.
  Defaults to FALSE as this is generally not of interest.}
}
\details{
  This function has been specifically written to read the grayscale tiff
  images which are produced by the Illumina scanners.  It is not generic
  enough to read all tiff files, although this functionality may be added
  in the future.
}
\value{
Returns a matrix with the same dimensions as the number of pixels in the
tiff file to be read in.
}

\author{Mike Smith}

\keyword{misc}

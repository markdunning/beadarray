\name{readTIFF}
\alias{readTIFF}
\title{Read the Illumina tiff images}
\description{
  Reads Illumina tiff images and produces a matrix of pixel values.
}
\usage{
  readTIFF(fileName, path = NULL, verbose = FALSE, xlim = NULL, ylim = NULL)
}

\arguments{
  \item{fileName}{String specifying the name of the tiff image to be read.}
  \item{path}{String specifying the path to the desired image.  The default value of NULL means the current working directory will be used.}
  \item{verbose}{If TRUE then details from the header of the tiff are
  printed as it is read.  These include things like the byte order, the
  number of pixels in the image, the number of tags in the header etc.
  Defaults to FALSE as this is generally not of interest.}
  \item{xlim}{Used to specify a subsection of the image to read in.  Takes a two entry vector e.g c(m,n) specifying that only pixels with x-coordinates between m and n (inclusive) should be read in.}
  \item{ylim}{Same a xlim, but for the y-cordinates.}
}
\details{
  This function has been specifically written to read the grayscale tiff
  images which are produced by the Illumina scanners.  It is not generic
  enough to read all tiff files, although this functionality may be added
  in the future.

  Given that the raw images can be quite large, functionality has also
  been included to read tiffs that have been compressed as either .bz2
  or .gz files.  Identification is performed based on the file extension
  and it is assumed that each tif is compressed individually.  Support
  for zip files may be added in the future.
}
\value{
Returns a matrix with the same dimensions as the pixels in the
tiff file to be read in.
}

\author{Mike Smith}

\keyword{IO}

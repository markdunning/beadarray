\name{BeadLevelNormalise}
\alias{BeadLevelNormalise}
\title{Bead Level Normalisation}
\description{
  A single function which combines the steps of reading bead level
  data, normalisation and summarising.
}
\usage{
BeadLevelNormalise(targets, imageManipulation="sharpen",
             backgroundCorrect="none", probes=NULL,imagesPerArray=2)
}
\arguments{
  \item{targets}{A targets object defining which bead level files to read} 
  \item{imageManipulation}{string specifying what image processing steps
    to use to get foreground intensities}
  \item{backgroundCorrect}{string specifying what background correction
    (if any) to do}
  \item{probes}{an optional list of ProbeIDs to use summarise}
  \item{imagesPerArray}{integer to specify how many images comprise each
    array. Used when we have to strips per array on the Human-Ref6 arrays}
}
\details{
  This is a memory and space efficient method for reading and normalising a
  large number of arrays at the bead level and based on the work of
  Henrik Bengsson. It comprises of a pre-processing, normalisation and
  probe summary.
  
  Quality control should be performed prior to this function is order to
  identify any outlier arrays and decide upon an image processing strategy.

  In the first stage of the function, we read each array in turn using
  the image processing steps specified by the imageManipulation and
  backgroundCorrect parameters and keep a running total of the value of
  each quantile of the distribution.

  For Human-6 arrays where each array is in fact two images, we read
  each image separately and bind together the resulting intensities
  before computing quantiles.

  Each array is then read a second time, normalised to the target
  distribution and then summarised.
}
\value{
ExpressionSetIllumina containing summarised data
}
\author{Mark Dunning}
\keyword{methods}
\name{beadRegistrationData-class}
\docType{class}
\alias{beadRegistrationData-class}
\alias{boxplot,beadRegistrationData-method}

\title{Class "beadRegistrationData" }

\description{
A class for storing information relating to the registration of the image.
}

\section{Slots/List Components}{
  Objects of this class contain the following slots
  \tabular{ll}{
  \code{layout}:\tab A list entry containing details of the structure of a BeadChip.  By default entries for the number of sections, the number of segments per section and a logical value stating whether the chip is two colour are created.  \cr
  \code{registrationData}:\tab A list with length equal to the total number of segments registered.  Each entry contains a vector of the differences between the within bead-type variances calculated using randomly assigned bead IDs and those calculated using the provided bead IDs.\cr
  \code{coordinateData}:\tab A list with length equal to the total number of segments registered.  Each entry contains the coordinates of the beads assigned to the segment, taken from the .locs file.\cr
  \code{cornerData}:\tab A list with length equal to the total number of segments registered.  Each entry contains a set of four coordinate pairs that mark the corners of the segment within it's image.  These can be extracted from the \code{coordinateData} and this slot may be removed at a later date. \cr
  \code{p95}:\tab A numeric vector storing the 95th percentile of bead intensities within each segment.\cr
  \code{imageLocations}:\tab Character vector storing the location of the tiff images of the array.  These are obtained from the sectionData slot of the \code{\link[=beadLevelData-class]{beadLevelData}} object that was processed by \code{\link{checkRegistration}}, and are required for visualising the bead locations. \cr
  \code{metrics}:\tab A data.frame, where each line contains the appropriate entry from the Metrics.txt file.  Currently not used, but may be later. \cr

  }
}

\section{Methods}{
  \describe{
    \item{boxplot(regScores, plotP95 = FALSE)}{}
  }
}

\seealso{\code{\link{checkRegistration}}}

\author{Mike Smith}
\keyword{classes}

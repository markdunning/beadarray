\name{displayTIFFImage}
\alias{displayTIFFImage}
\alias{showTIFF}
\alias{FindIlluminaSpotl}
\alias{Recenter}
\alias{IdentifyIntensity}
\alias{IdentifyBeads}

\title{ Display part of a TIFF Image}
\description{
Function for looking at a particular part of an image scanned from a BeadArray
}

\usage{
displayTIFFImage(BLData, array, a = 680:720, b = 680:720, flip=TRUE,
showOutliers=TRUE, locateBeads=FALSE, showUnregistered=FALSE, outliers=NULL)
}

\arguments{
  \item{BLData}{an BeadLevelList object containing bead level data}
  \item{array}{numeric value for the array that we want to look at}
  \item{a}{vector giving the range of x values for the region we want to look at}
  \item{b}{vector giving the range of y values for the region we want to look at }
  \item{flip}{if TRUE then the input y range of coordinates will be flipped so that y=0 will be at the top of the array}
  \item{showOutliers}{if TRUE then the location of outliers within the region will be highlighted on the plot}
  \item{locateBeads}{if TRUE then the plot will be interactive}
  \item{showUnregistered}{if TRUE then any unregistered points inside the region will be highlighted} 
 \item{outliers}{List of BeadIDs which have been discovered to be
   outliers previously}
}
\details{
Loads up the original pgm image scanned from the array and displays the part of the image defined 
by a and b. Intensities of pixels are given on the log2 scale with a brighter green indicating a value near 16.

The 'showOutliers' parameter sets if outliers in the region are highlighted in the region with a default value of TRUE. All outliers which have intensity higher or lower than the mean for their bead type will 
be displayed as blue or red dots respectively. Any beads which have negative intensity will be displayed as yellow dots.

The 'showUnregistered' parameter allows beads which were not correctly decoded to be highlighted. Setting this parameter to TRUE
causes any unregistered beads in the region to be displayed as grey dots.

The 'flip' parameter arises from the fact that when we plot in R, we have the line y=0 at the \textit{bottom} of the plot, whereas on BeadArrays the y values are counted from the \textit{top} of the array. 

Setting the locateBeads parameter to TRUE will make the plots produced interactive. In this mode, clicking on outliers on the image will display information about the unsharpened intensity, 
the local background level along with the probeID for the bead and mean of all beads on the array with that probeID.

}
\value{
Plot is produced on the current graphical device.
}

\author{Mark Dunning and Isabelle Camilier}

\seealso{
\code{\link{getUnregisteredBeads}}
}

\examples{

data(BLData)

##The image files that are to be read need to be in the working directory

##displayTIFFImage(BLData, 1, a=1000:1400, b=1200:1400)

##To display unregistered beads

##displayTIFFImage(BLData, 1, a=1300:1600, b=400:600, showUnregistered=TRUE)

##To use plot interactively and find information about beads clicked on set 'locateBeads' parameter

##displayTIFFImage(BLData, 1, a=1300:1600, b=400:600, locateBeads=TRUE)

}
\keyword{hplot}

\name{ArrayMask}
\alias{addArrayMask}
\alias{removeArrayMask}
\alias{clearArrayMask}
\alias{showArrayMask}

\title{Array Mask}
\description{
Functions to edit or display array masks.
}

\usage{
addArrayMask(BLData, array, SAM = FALSE, nrow = 50, ncol = 50, high = "red",
	     low = "yellow", zlim = c(7,15), override = FALSE)
removeArrayMask(BLData, array, SAM = FALSE, nrow = 50, ncol = 50,
                high = "red", low = "yellow",
                zlim = c(7,15), override = FALSE)
showArrayMask(BLData,array, SAM = FALSE, elim = TRUE, override = FALSE)
clearArrayMask(BLData,array)
}

\arguments{
  \item{BLData}{A BeadLevelList object.}
  \item{array}{The number of an array in the BeadLevelList object.}
  \item{SAM}{Logical. If \code{TRUE}, display a hexagonal overlay where appropriate.}
  \item{elim}{Logical. If TRUE, plot eliminated beads (with blue crosses).}
  \item{nrow, ncol, high, low, zlim}{Arguments passed to \code{imageplot} - see the help for \code{imageplot} for details.}
  \item{override}{Logical. Plotting a large mask can cause slowdown problems. By default, if more than 200 000 beads are masked, the current mask will not be plotted. You can force the mask to be plotted by setting this argument to \code{TRUE}, however beware as this may cause slower systems to freeze.}
}

\details{
These functions are used to manipulate the mask on a single array.

\code{addArrayMask} adds beads in a specified region to the mask (i.e. sets their weights to 0).

\code{removeArrayMask} removes beads in a specified region from the mask
(i.e. sets their weights to 0)
 

On calling either of these functions, an imageplot is displayed. Click
on the plot to define vertices of a polygon, in order. Having specified
the last vertex, right-click to close the polygon. A plot is then
produced of the beads currently masked on the array (in grey) and the
polygon just defined (in red), with a menu prompting you to accept the
displayed region - if you do so, then all beads in the polygon then be
masked or unmasked as appropriate. Alternatively, you can right-click on
the image without defining any vertices, thus leaving the mask
unchanged. 
 

If making the change would result in all beads of a certain probe ID
being completely covered by the mask, then the functions return a
warning message, and the beads eliminated in this way are highlighted
with blue crosses on the plot.
 

\code{showArrayMask} plots the beads on an array which have been masked,
over a plot of outliers.
 

\code{clearArrayMask} clears the mask on an array, removing all weights
associated with it.
 
}

\value{
  None returned
}


\author{Jonathan Cairns}

\examples{
##data(BLData)
##addArrayMask(BLData,1)
##showArrayMask(BLData,1)
}

\keyword{misc}

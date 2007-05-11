\name{readIllumina}
\alias{readIllumina}

\title{Read BeadArray information into R}
\description{
Uses csv or txt and image files (where available) to load information
about each array used in a BeadChip or SAM experiment.
}
\usage{
readIllumina(arrayNames=NULL, path=".", textType=".csv", annoFile=NULL,
        targets=NULL, imageManipulation="sharpen", backgroundSize=17,
        storeXY=TRUE, sepchar="_", metrics=FALSE,
        metricsFile="Metrics.txt", backgroundMethod="none", offset=0,
        normalizeMethod="none", ...)
}

\arguments{
  \item{arrayNames}{character vector containing names of arrays to be
    read in (if present)}
  \item{path}{character string specifying the location of files to be
    read by the function}
  \item{textType}{character string specifiying the extension of the
    files which store the xy information (typically \code{".txt"} or
    \code{".csv"})}
  \item{annoFile}{the name of a file (such as \code{.opa}) containing
    annotation information about each bead type}
  \item{targets}{data.frame containing sample information}
  \item{imageManipulation}{Currently the only options here are "none"
    or "sharpen", which indicates the Illumina sharpening mask will be used
    prior to foreground averages being calculated}
  \item{backgroundSize}{numeric value which defines the size of the n x
    n box we will use to calculate local background values}
  \item{storeXY}{logical scalar, indicating whether the xy coordinates
    should be stored}
  \item{sepchar}{character string which separates row and column in file
    name (default value is \code{"_"})}
  \item{metrics}{logical scalar, indicating whether the scanner metrics
    file \code{metricsFile} is to be read in}
  \item{metricsFile}{name of the scanner metrics file}
  \item{backgroundMethod}{method to use for background correcting the
    data.  Options are \code{"none"}, \code{"subtract"}, \code{"half"},
    \code{"minimum"}, \code{"edwards"}, \code{"normexp"} or \code{"rma"}}
  \item{offset}{numeric value to add to intensities}
  \item{normalizeMethod}{method to use to normalize the background
    corrected bead-level data. Options are \code{"none"},
    \code{"quantile"} and \code{"vsn"}. Note that the normalization occurs
    at the bead-level and is only available for two-colour data at this
    stage}
  \item{...}{other arguments}
}

\details{
This function can be used to read in bead level information using the
raw ./tif and csv ot txt files output by the scanner.
Where sample information is available, the \code{targets} data.frame can
be read in using 'read.table'.

The \code{annoFile} contains information on each bead type on the array.

The csv or txt files are used to specify the location and identity of
each bead on the array and must contain columns for the $x$ and $y$
postion of each bead as well as a ProbeID. For two-colour arrays, this
information is required for each channel.

The foreground and background intensities of each bead are then
calculated and stored separately. For foreground calculations the
sharpening mask used by Illumina is used prior to an averaging over
the 9 pixels in a 3 x 3 square closest to the bead centre (we can choose
to calculate the values obtained without this mask). For the background
intensity we look at the pixels in a 17 x 17 square and calculate the
mean of the 5 lowest values inside this square.  If a bead is
particularly close to the edge of the image then it isn't possible to
take a 17 x 17 square for sampling.  In this case the background value
is taken as the mean of all the previous background values and the user is
informed that this has been done.

The function creates an BeadLevelList containing foreground and
background intensities for each bead on each array.
}
\value{
  BeadLevelList object
}

\author{Mark Dunning, Mike Smith}

\examples{
#BLData = readIllumina()
 
#targets = read.table("targets.txt", header=T)

#targets
#May take a while to run
#BLData.s = readIllumina(arrayNames=target$Institute.Sample.Label, targets=targets)

#Create foreground intensities without using sharpening. Should take less time
#BLData.ns = readIllumina(arrayNames=targets$Institute.Sample.Label, targets=targets, sharpen=FALSE)

#Use a 9 x 9 are around each bead to calculate the background
#BLData.ns.9 = readBeadImages(arrayNames=target$Institute.Sample.Label, targets=targets, sharpen=FALSE, backgroundSize=9)
}
\keyword{IO}

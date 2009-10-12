\name{readIllumina}
\alias{readIllumina}

\title{Read bead-level data into R}
\description{
Uses .csv or .txt and TIFFs (where available) to load information
about each bead on each array in a BeadChip or SAM experiment.
}
\usage{
readIllumina(arrayNames=NULL, path=".", textType=".txt",
        annoPkg=NULL, useImages=TRUE,  
        singleChannel=TRUE, targets=NULL, 
        imageManipulation="sharpen", backgroundSize=17, backgroundCalc = "mean",
        storeXY=TRUE, sepchar="_", dec=".", metrics=FALSE,
        metricsFile="Metrics.txt", backgroundMethod="subtract", 
        offset=0, normalizeMethod="none",
        tiffExtGrn="_Grn.tif", tiffExtRed="_Red.tif", ...)
}

\arguments{
  \item{arrayNames}{character vector containing names of arrays to be
    read in.  If \code{NULL}, all arrays that can be found in the current 
    working directory will be read in.}
  \item{path}{character string specifying the location of files to be
    read by the function}
  \item{textType}{character string specifiying the extension of the
    files which store the bead-level information.  Typically \code{".txt"}, 
    \code{".csv"} or \code{"_perBeadFile.txt"}.}
  \item{annoPkg}{character string specifying the annotation package 
     for the arrays being read in (only available for certain 
     expression arrays at present).  Default value is
     \code{"illuminaProbeIDs"} which is not an annotation package, and
     indicates that Illumina bead IDs have been used to identify each bead.}
  \item{useImages}{logical. If TRUE, the foregound and background 
    values are retrieved from the TIFFs.  When FALSE, the 
    intensity values in the text files are used.  Note that
    background values will not be available (set to 0) when FALSE, 
    as the current option in BeadScan is to store background
    corrected intensities.}
  \item{singleChannel}{logical.  Set to TRUE if the data is single 
   channel (Green images only) or FALSE for two-colour (both Green and Red 
   data available).}
  \item{targets}{\code{data.frame} containing sample information}
  \item{imageManipulation}{character string specifying which image analysis 
    method to use.  The current options are \code{"none"} (no image manipulation 
    or \code{"sharpen"} (the Illumina sharpening mask will be used prior to 
    the foreground averages being calculated).}
  \item{backgroundSize}{integer value which defines the size of the n x
    n region (in pixels) used to calculate local background values.  
    Default value is 17}
  \item{backgroundCalc}{Specifies how the background value for each bead is calculated.
    Currently the only options are "mean" and "median", which return the specified statistic 
    for the five dimest pixels in the region defined by the \code{backgroundSize} argument.  
    If any other value is given then the mean is returned.}
  \item{storeXY}{logical scalar, indicating whether the xy coordinates
    should be stored}
  \item{sepchar}{character string which separates the row and column positions
    in the file names (default value is \code{"_"})}
  \item{dec}{character used in the files for decimal points.  The default
    value is \code{"."}}
  \item{metrics}{logical scalar, indicating whether the scanner metrics
    file \code{metricsFile} is to be read in}
  \item{metricsFile}{name of the scanner metrics file (\code{"Metrics.txt"}
    by default)}
  \item{backgroundMethod}{method to use for background correcting the
    data.  Options are \code{"none"}, \code{"subtract"}, \code{"half"},
    \code{"minimum"}, \code{"edwards"}, \code{"normexp"} or \code{"rma"}}
  \item{offset}{numeric value to add to intensities}
  \item{normalizeMethod}{method to use to normalize the background
    corrected bead-level data. Options are \code{"none"},
    \code{"quantile"} and \code{"vsn"}. Note that the normalization occurs
    at the bead-level and is only available for two-colour data at this
    stage}
  \item{tiffExtGrn}{character string specifying the file extension of the Cy3 (Green)
    images.  Default is \code{"_Grn.tif"}}
  \item{tiffExtRed}{character string specifying the file extension of the Cy5 (Red)
    images (where present).  Default is \code{"_Red.tif"}}
  \item{...}{other arguments}
}

\details{
This function can be used to read in bead-level information from the
raw .tif and .csv or .txt files output by BeadScan.

Note that the .txt or .csv files must contain the raw data for each bead 
on each array/strip, not the summarised data.
The .csv or .txt files specify the location and identity of each bead on 
the array and must contain columns for the $x$ and $y$ position of each 
bead as well as a ProbeID.  For two-colour arrays, this information is 
required for each channel.

The foreground and background intensities of each bead are 
calculated from the images when \code{useImages=TRUE}.  For the 
foreground calculations the sharpening mask used by Illumina is 
applied prior to averaging over the 9 pixels in a 3 x 3 square closest 
to the bead centre by default (\code{imageManipulation="sharpen"}). 
If \code{imageManipulation="none"}, no sharpening is performed. 
The local background estimate for each bead is calculated by averaging 
the 5 minimum pixels in a 17 x 17 square around each bead centre from the 
unsharpened image.  The size of this window is controlled by the 
\code{backgroundSize} argument.  If a bead is too close to the edge of 
the image, it is ignored.

When \code{useImages=FALSE}, the intensities from the .txt or .csv files are
stored as the foreground values for each bead.  Note that the values stored in 
these files have already undergone a local background adjustment, so the 
background values are set to 0.

To keep track of the samples hybridised to each array, we recommend 
using a \code{targets} \code{data.frame}, which lists each strip/array in the rows,
and experimental information about each strip/array in the columns 
(sample, array name, etc.)
Targets information can easily be created and saved in tab delimited text format,
read in using \code{read.table} and passed to \code{readIllumina} using the  
\code{targets} argument.

%Note that for BeadChips which have multiple strips per array, it is necessary 
%to have a row in the targets file for each strip, as there will be a .tif and 
%.txt or .csv file for each strip.  

The pairs of strips from a BeadChip can be combined when the data is summarised with \code{createBeadSummaryData}.

The function creates a \code{BeadLevelList} containing foreground and
background intensities for each bead on each array.

NOTE:  Reading in bead-level data, particularly with the TIFFs is memory intensive.  
For example, reading in text and image data from a Human-6 BeadChip requires several 
Gigabytes of RAM.  If you have limited memory, it is recommnded that you read in the data using the \code{useImages=FALSE} option.
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
#BLData.s = readIllumina(arrayNames=target$Institute.Sample.Label, targets=targets, imageManipulation="none")

#Create foreground intensities without using sharpening. Should take less time
#BLData.ns = readIllumina(arrayNames=targets$Institute.Sample.Label, targets=targets, imageManipulation="sharpen")

}
\keyword{IO}

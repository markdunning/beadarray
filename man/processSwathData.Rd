\name{processSwathData}
\alias{processSwathData}
\title{Prepare iScan data for use with beadarray}
\description{
Data from Illumina's newer iScan system come in a different format to the previous BeadScan data.  This function is intended to transform the data into a format compatible with beadarray.
}


\usage{
processSwathData(inputDir = NULL, outputDir = NULL, twoColour=NULL, textstring="_perBeadFile.txt", segmentHeight = 326, segmentWidth = 397, fullOutput = FALSE, newTextString=".txt", verbose = FALSE)
}

\arguments{
\item{inputDir}{Character string specifying the directory containing the data to be processed.  If left \code{NULL} this defaults to the current working directory.}
\item{outputDir}{Similar to the above, specifying the directory where the output should be written.  If left \code{NULL} this defaults to the current working directory.}
\item{twoColour}{Boolean value specifying whether the data is one or two channel.  If left NULL the function will attempt to determing the number of channels be examing the files present in the directory.}
\item{textstring}{String specifying the suffix to identify the text file containing original bead-level data.}
\item{segmentHeight}{Each array section is made up of several segments of beads, arranged in an hexagonal pattern.  This value specifies the number of rows in a segment, and can be found in the accompanying .sdf file (<SizeGridX>).  In the future we will attempt to automate this.}
\item{segmentWidth}{Similar to the above argument, this specifies the number of columns in an array segment (<SizeGridY>).}
\item{fullOutput}{Boolean value specifying the type of output.  More details are given below.}
\item{newTextString}{Suffix for the two new "per swath" bead-level text file.}
\item{verbose}{Boolean value that, if TRUE, directs the function to print progress to the screen.}
}

\details{
Data from the iScan system comes with two images of each array section (along with two .locs files), which are labelled Swath1 and Swath2.  These two images are of the two halves of the array section, with an overlapping region in the middle.  However, there is only one bead-level text file, with no indication as to which of the two images each entry comes from.  Given this, simply reading the bead-level text file will result in any function that uses bead locations performing undesirably.

This function works to try and deconvolute the bead-level data and create two files, one per swath, which can then be read independently into beadarray.

The exact content of the output files depends upon the \code{fullOutput} argument.  If the default value of \code{FALSE} is selected the function compares the coordinates in the "perBeadFile.txt" file with those found in the two .locs files, in order to determine which swath each bead is from.  From this two new bead-level text files (with names containing "Swath1" and "Swath2" are created, containing all the beads from the orginal file.  If the .locs files are not present this process will fail and the assigning of beads to swaths cannot be performed.  

If the \code{fullOutput} argument is set to \code{TRUE} both text files contain all the beads that can be identified in their respective images.  Assuming that both the .locs and .tif files are present, we can use the bead-centre coordinates stored in the .locs files to calculate intensity values for beads in both images, even if no intensity is recorded Illumina's text file.  Any bead that is found in the overlapping region, and thus appearing in both images, will have two intensities calculated.  The files that are created contain the same data as regular bead-level text files, but with an additional column entitled "Weights".  Beads that appear in both images are assigned a weight of 0.5 (since the are two intensity values for them), whilst all other beads are given a weigth of 1.  If the TIFF images aren't available the intensities cannot be calculate, so the output will default to the same as if \code{fullOutput = FALSE}.
}
\value{
This function is called for its side effects, which is to produce two text files containing the beads match to the approriate swath.  No value it returned by the function.
}
\author{
Mike Smith and Andy Lynch
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{IO}


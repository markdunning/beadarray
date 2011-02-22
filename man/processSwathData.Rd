\name{processSwathData}
\alias{processSwathData}
\title{Prepare iScan data for use with beadarray}
\description{
Data from Illumina's newer iScan system come in a different format to the previous BeadScan data.  This function is intended to transform the data into a format compatible with beadarray.
}


\usage{
processSwathData(twoColour=TRUE, textstring="_perBeadFile.txt", segmentHeight = 326, segmentWidth = 397, swathOverlap = 81, fullOutput = TRUE, newTextString=".txt", verbose = TRUE)
}

\arguments{
\item{twoColour}{Boolean value specifying whether the data is one or two channel.}
\item{textstring}{String specifying the suffix to identify the text file containing original bead-level data.}
\item{segmentHeight}{Each array section is made up of several segments of beads, arranged in an hexagonal pattern.  This value specifies the number of rows in a segment, and can be found in the accompanying .sdf file.  In the future we will attempt to automate this.}
\item{segmentWidth}{Similar to the above argument, this specifies the number of columns in an array segment.}
\item{swathOverlap}{Specifys the number of columns that appear in both swaths.  Again this will hopefully be automated in the future.}
\item{fullOutput}{Boolean value specifying the type of output.}
\item{newTextString}{Suffix for the two new "per swath" bead-level text file.}
\item{verbose}{Boolean value that, if TRUE, directs the function to print progress to the screen.}
}

\details{
Data from the iScan system comes with two images of each array section (along with two .locs files), which are labelled Swath1 and Swath2.  These two images are of the two halve of the array section, with an overlapping region in the middle.  However, there is only one bead-level text file, with no indication as to which of the two images each entry comes from.  Given this, simply reading the bead-level text file will result in any function that uses bead locations to perform undesirably.

This function works to try and deconvolute the bead-level data and create two files, one per swath, which can then be read independently into beadarray.
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


\name{processSwathData}
\alias{processSwathData}
\title{Prepare iScan data for use with beadarray}
\description{
Data from Illumina's newer iScan system come in a different format to the previous BeadScan data.  This function is intended to transform the data into a format compatible with beadarray.
}


\usage{
processSwathData(twocolour=TRUE, textstring="_perBeadFile.txt", section_height = 326, section_width = 397, swath_overlap = 81, fullOutput = TRUE, newTextString=".txt", verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{twoColour}{Boolean value specifying whether the data is one or two channel.}
\item{textstring}{String specifying the suffix to identify the text file containing original bead-level data.}
\item{section_height}{}
\item{section_width}{}
\item{swath_overlap}{}
\item{fullOutput}{}
\item{newTextString}{Suffix for the two new "per swath" bead-level text file.}
\item{verbose}{Boolean value that, if TRUE, directs the function to print progress to the screen.}
}

}
\details{
Data from the iScan system comes with two images of each array section (along with two .locs files), which are labelled Swath1 and Swath2.  These two images are of the two halve of the array section, with an overlapping region in the middle.  However, there is only one bead-level text file, with no indication as to which of the two images each entry comes from.  Given this, simply reading the bead-level text file will result in any function that uses bead locations to perform undesirably.

This function works to try and deconvolute the bead-level data and create two files, one per swath, which can then be read independently into beadarray.
}
\value{

}
\author{
Mike Smith and Andy Lynch
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{IO}


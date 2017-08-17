\name{readLocsFile}
\alias{readLocsFile}
\title{
    Read ``.locs'' file.
}
\description{
    Reads the binary Illumina bead location files and returns a matrix of
    the coordinate pairs for every bead on the array.
}
\usage{
    readLocsFile(fileName)
}

\arguments{
    \item{fileName}{
        A string containing the name of the ``.locs'' file to be read.
    }
}

\details{
  The locs file contains bead centre locations for every bead on the
  array, unlike the bead level text files, with contain just the beads
  that were decoded. Reading these can be useful if one wants to verify
  that the image registration was successful, or is interested in the
  locations of the undecoded beads.

  The locs file itself is in a binary format, with each of the bead
  locations stored as a pair of doubles.  The first 2 bytes contain
  header information, with the 3rd byte containing the number of probes
  on the array.  The location information begins with the 4th byte.
}

\value{
    Returns a two column matrix of bead coordinates, one row per bead.
}

\author{
    Mike Smith
}

\keyword{IO}


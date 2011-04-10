\name{checkRegistration}
\alias{checkRegistration}
\title{
Perform check for misregistered array segments.
}
\description{
Occasionally the registration of an array can go wrong, with the bead centres found in the wrong place in an image.  The effective result of this is a scrambling of the bead IDs.  In order to check for this we can examine the within bead-type variance across the array.  In cases where registration has failed we would expect to see a large jump in this value compared to correctly registered arrays.  This function computes this statistic for each array segment (since each segment is registered independently) and returns them to the user for inspection.

Note that the function requires that the sdf file and locs file are present, and has particular expectations towards their file names and locations.
}
\usage{
    checkRegistration(BLData, array = 1)
}

\arguments{
    \item{BLData}{An object of class \code{\link{beadLevelData-class}}.}
    \item{array}{Integer specifying the index of the arrays to be checked.  Can be a vector to process multiple arrays e.g. \code{1:12}.}
}

\value{
Returns a list with an entry for each section specified by the \code{array} argument.  Within each of these entries contains a vector storing the mean within bead-type for each segment of the array.  If the array is two-colour then there will be seperate entries for the green and red channel as these are registered separately by the scanner.
}

\author{
Mike Smith
}

\references{
Smith ML, Dunning MJ, Tavare S, Lynch AG. Identification and correction of previously unreported spatial phenomena using raw Illumina BeadArray data. BMC Bioinformatics (2010) 11:208   
}

\keyword{misc}


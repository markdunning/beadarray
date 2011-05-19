\name{checkRegistration}
\alias{checkRegistration}
\title{
Perform check for misregistered array segments.
}
\description{
Occasionally the registration of an array can go wrong, with the bead centres found in the wrong place in an image.  The effective result of this is a scrambling of the bead IDs.  

Note that the function requires that the sdf file and locs file are present, and has particular expectations towards their file names and locations.
}
\usage{
    checkRegistration(BLData, array = 1)
}

\arguments{
    \item{BLData}{An object of class \code{\link{beadLevelData-class}}.}
    \item{array}{Integer specifying the index of the arrays to be checked.  Can be a vector to process multiple arrays e.g. \code{1:12}.}
}

\details{
    In order to check for mis-registration we can examine the within bead-type variance across the array.  This function computes this statistic twice for each array segment (since each segment is registered independently), once using the given bead IDs and once using a randomly assigned set of IDs.  The former is then subtracted from the later.  In cases where the registration has worked successfully we expect the majority of these values to be greater than zero, which for misregistered arrays the differences should be centred about zero.
}

\value{
Returns an object of class \code{\link[=beadRegistrationData-class]{beadRegistrationData}}.
}

\author{
Mike Smith
}

\references{
Smith ML, Dunning MJ, Tavare S, Lynch AG. Identification and correction of previously unreported spatial phenomena using raw Illumina BeadArray data. BMC Bioinformatics (2010) 11:208   
}

\keyword{misc}


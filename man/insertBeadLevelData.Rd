\name{insertBeadLevelData}
\alias{insertBeadLevelData}
\title{
    Add or modify data in a beadLevelData object
}
\description{
    Add or modify data in a \code{\link{beadLevelData}} object.
}
\usage{
    insertBeadLevelData(BLData, array = 1, what, data)
}

\arguments{
    \item{BLData}{
        An object of class \code{\link{beadLevelData}}.
    }
    \item{array}{
        Positive integer specifying what section should be modified.
    }
    \item{what}{
        Name of the data that is being modified.  It `what' doesn't exist then a new entry is created using the name specified in this argument.
    }
    \item{data}{
        A numeric vector to be stored, the same length as the number of beads in the section specified by the array argument.
    }
}

\details{
This function allows the \code{beadData} slot of the \code{\link{beadLevelData}} object to be modified for a given array. 

}

\value{
    Returns an object of class \code{\link{beadLevelData}}.
}

\author{
    Mike Smith
}

\keyword{IO}


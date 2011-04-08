\name{insertBeadData}
\alias{insertBeadData}
\alias{removeBeadData}
\title{
    Add, modify or remove data in a beadLevelData object
}
\description{
    Add, modify or remove data in a \code{\link[=beadLevelData-class]{beadLevelData}} object.
}
\usage{
    insertBeadData(BLData, array = 1, what, data)
    removeBeadData(BLData, array = 1, what)
}

\arguments{
    \item{BLData}{
        An object of class \code{\link{beadLevelData-class}}.
    }
    \item{array}{
        Positive integer specifying what section should be modified.
    }
    \item{what}{
        Name of the data that is being modified.  If `what' doesn't exist then a new entry is created using the name specified in this argument.
    }
    \item{data}{
        A numeric vector to be stored, the same length as the number of beads in the section specified by the array argument.
    }
}

\details{
These functions allow the \code{beadData} slot of the \code{\link{beadLevelData-class}} object to be modified for a given array. 
}

\value{
    Returns an object of class \code{\link{beadLevelData-class}}.
}

\author{
    Mike Smith
}

\keyword{IO}


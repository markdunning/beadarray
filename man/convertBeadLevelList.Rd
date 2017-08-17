\name{convertBeadLevelList}
\alias{convertBeadLevelList}
\title{
    Convert a BeadLevelList object into a beadLevelData object
}
\description{
    As of beadarray version 2.0 the BeadLevelList class has been deprecated and replaced by the beadLevelData class.  Whilst these are superficially similar, the way the data are stored is quite different, meaning most functionality within the package is no longer compatible with the original BeadLevelList class.
    
    This function converts any object that is of the old BeadLevelList class into a beadLevelData object.
}
\usage{
    convertBeadLevelList(BeadLevelList)
}

\arguments{
    \item{BeadLevelList}{
        An object of class BeadLevelList
    }
}

\value{
    Returns an object of class beadLevelData.
}

\author{
    Mike Smith
}

\seealso{
    \code{\link{beadLevelData-class}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{IO}


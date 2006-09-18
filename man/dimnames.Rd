\name{dimnames}
\alias{dimnames.BeadLevelList}

\title{Retrieve the Dimension Names of an BeadLevelList or BeadSummaryList Object}
\description{
Retrieve the dimension names of a microarray data object.
}
\usage{
\method{dimnames}{BeadLevelList}(x)
}
\arguments{
  \item{x}{an object of class \code{BeadLevelList}}
}
\details{
The dimensionnames of an microarray object are the same as those of the most important matrix component of that object.

A consequece is that \code{rownames} and \code{colnames} will work as expected.
}
\value{
Either \code{NULL} or a list of length 2.
If a list, its components are either 'NULL' or a character vector the length of the appropriate dimension of \code{x}.
}
\author{Gordon Smyth, edited by Mike Smith}
\seealso{
  \code{\link[base]{dimnames}} in the base package.
}
\keyword{array}

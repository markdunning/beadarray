\name{setWeights}
\alias{setWeights}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set weights from BASH
}
\description{
Function for committing the weights calculated by BASH into the beadLevelData object
}
\usage{
setWeights(BLData, wts, array, combine = FALSE, wtName = "wts")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{
a \code{beadLevelData} object
}
  \item{wts}{
the wts component of the BASH output
}
  \item{array}{
a vector of arrays that we want to set the weights for
}
  \item{combine}{
if TRUE combine the weights with existing weights (if they exist)
}
  \item{wtName}{
name of column to assign weights to
}
}

\value{
Modified beadLevelData object

}
\author{
Mark Dunning
}



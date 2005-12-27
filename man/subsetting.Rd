\name{subsetting}
\alias{[.BeadLevelList}
\alias{[.BeadSummaryList}
\alias{printHead}
\title{Subset BeadLevelList, BeadSummaryList}
\description{
Extract a subset of an \code{BeadLevelList} or \code{BeadSummaryList}
}
\usage{
object[i, j]
}
\arguments{
  \item{object}{object of class \code{BeadLevelList} or \code{BeadSummaryList}}
  \item{i,j}{elements to extract. \code{i} subsets the genes or spots while \code{j} subsets the arrays}
}
\details{
\code{i,j} may take any values acceptable for the matrix components of \code{object}.
See the \link[base]{Extract} help entry for more details on subsetting matrices.
}
\value{
An object of the same class as \code{object} holding data from the specified subset of genes and arrays.
}
\author{Gordon Smyth, modified by Mike Smith for BeadSummaryList object}
\seealso{
  \code{\link[base]{Extract}} in the base package.
%  \link{3.ReadingData} gives an overview of data input and manipulation functions in LIMMA.
}

\keyword{manip}

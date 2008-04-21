\name{rankIvariantNormalise}
\alias{rankInvariantNormalise}
\title{Rank Invariant normalise data in a matrix}
\description{
Normalise expression matrix using rank invariant genes
}
\usage{
rankInvariantNormalise(exprs, T)
}
\arguments{
  \item{exprs}{a matrix of expression values}
  \item{T}{A target distribution vector to normalise the data to.
    The default is NULL in which case the average quantiles are used}
}
\details{
  Using the \code{normalize.invariantset} function from the affy package, 
  we find a list of rank invariant genes whose rank does not change significantly 
  between the columns of \code{exprs}.  We then fit a normalising curve to each 
  array using the values of the rank invariant genes of the array and a target
  distribution.

  The target distribution may be specified by the user and by default
  is the vector of average quantiles across all arrays.
}

\value{
Matrix of normalised expression data with the same dimensions as \code{exprs}.
}
\author{Mark Dunning}

\examples{
data(BSData)

BSData.ri  = assayDataElementReplace(BSData, "exprs", rankInvariantNormalise(exprs(BSData)))
}
\keyword{methods}
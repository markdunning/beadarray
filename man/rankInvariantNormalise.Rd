\name{rankIvariantNormalise}
\alias{rankInvariantNormalise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Rank Invariant Normalisation for Bead Summary data}
\description{
Normalise expression matrix using rank invariant genes
}
\usage{
rankInvariantNormalise(exprs, T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exprs}{An expression matrix from a ExpressionSet object} 
  \item{T}{A target distribution vector used to normalise to. The default is NULL in which case the average quantiles are used}
}
\details{

Using the normalize.invariantset from affy, we find a list of rank invariant genes whose rank does not change significantly between the columns of exprs. We then 
fit a normalising curve to each array using the values of the rank invariant genes of the array and a target distribution.

The target distribution may be specified by the user and by default is the vector of average quantiles across all arrays.

}

\value{
Expression matrix with the same dimensions as the input

}
\author{Mark Dunning}


\examples{

data(BSData)

BSData.ri  = assayDataElementReplace(BSData, "exprs", rankInvariantNormalise(exprs(BSData)))



}

\keyword{methods}
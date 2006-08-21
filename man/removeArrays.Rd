\name{removeArrays}
\alias{removeArrays}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Finds unregistered beads}
\description{
Remove a set of arrays from a BeadSummaryList 
}
\usage{
removeArrays(BSData, array)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BSData}{an BeadSummaryList containing bead summary data}
  \item{array}{numeric vector giving the indices of the arrays to be removed}
}
\details{
Arrays can removed from an BeadSummaryList (bead summary data) object if they are found to contain bad data. Caution:- once arrays have been removed they cannot be recovered. It is recommended to create a new BeadSummaryList object for the output of this function.

}
\value{
  BeadSummaryList object
  
}

\author{Mark Dunning}

\examples{
data(BSData)
BSData$G[1,]
BSData.new=removeArrays(BSData, array=1:2)
BSData.new$G[1,]
} 
\keyword{methods}
\name{getVariance}
\alias{getVariance}
\title{Gets the bead-type variances from an ExpressionSetIllumina Object}
\description{
Calculates the variance for each bead-type on each array from an \code{ExpressionSetIllumina} object.
}
\usage{
getVariance(object, offset=0)
}
\arguments{
  \item{object}{\code{ExpressionSetIllumina} object}
  \item{offset}{numeric value to add to the variances to avoid very small values}
  }
\details{
  \code{getVariance} uses the \code{se.exprs} and \code{NoBeads} slots in \code{assayData} to calculate the variances for each bead-type on each array.
   }
\value{
A matrix containing the variances.
}
\author{Matt Ritchie}
\examples{
data(BSData)
v = getVariance(BSData)
boxplot(as.data.frame(log2(v)), ylab="log2var", xlab=colnames(BSData), las=2)
}
\keyword{manip}

\name{histBeadCounts}
\alias{histBeadCounts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot numbers of each bead}
\description{
 Function to plot on a histogram the numbers of each bead type found on an array.
}
\usage{
histBeadCounts(BLData, array, main = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{BeadLevelList object}
  \item{array}{numeric value giving the array we want to plot for}
  \item{main}{specify a title for the plot}  	
  \
}
\details{
 It is known that the number of each bead type to be found on an array follows a Poisson distribution with a mean of around 30. This plot helps us
to see how well this assumption is met. If too many bead types have less than 30 beads it might indicate a problem with the array and compromise the 
robustness of the experiment.

We show a Poisson distribution with mean ~50,000 / 1500 on top of the histogram so it can be seen how different the numbers of beads we observe is to
how many we expect to see.

WARNING: THIS FUNCTION WILL NOT BE EFFICIENT FOR WHOLE GENOME DATA

}
\value{
Histogram is plotted.
}

\author{Mark Dunning}

\examples{

data(BLData)

histBeadCounts(BLData, 1)

}
\keyword{hplot}
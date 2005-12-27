\name{getMeanIntensities}
\alias{getMeanIntensities}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find Bead Type Averages}
\description{
Function uses the bead summary data to find the average intensity of a particular bead type on each array in the experiement.
}
\usage{
getMeanIntensities(BSData, probeID)
}

\arguments{
  \item{BSData}{an BeadSummaryList object containing bead level data}
  \item{probeID}{numeric value giving the probeID we want values for}
}
\details{
The function returns the calculated bead summary value that has been calculated for each array in the experiment
}
\value{
vector of length n (where n is number of arrays in experiment) giving the average intensity of the bead type on each array in the experiment.
}

\author{Mark Dunning}

\seealso{\code{\link{createBeadSummaryData}}}

\examples{

data(BSData)

##To get the average intensitiy of bead with probeID 2 on all 4 arrays in the experiement use:

getMeanIntensities(BSData, 2)

}
\keyword{methods}
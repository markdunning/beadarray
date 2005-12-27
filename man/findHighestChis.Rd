\name{findHighestChis}
\alias{findHighestChis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find least randomly distributed beads}
\description{
 Function for finding which bead types on an array are least randomly positioned.

}
\usage{
findHighestChis(BLData, array, limit = 14)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{an BeadLevelList object containing bead level data}
  \item{array}{numeric value for which array we want to use}
  \item{limit}{numeric value which determines the threshold value for the chi-statistic.}
}
\details{
The $\chi^2$ statistic is defined as:

\begin{center}

$\chi^2 = \sum^k_{i=1} (O_i - E_i)^2 / E_i$
\end{center}

where O and E are vectors of length 8 defining the observed and expected number of beads in each section respectively. 
This provides a means of systematically testing the distribution of every bead type on a given array and we can record each bead type with a sufficiently high value of $\chi^2$.

If the value for a particular bead type is found
to be greater than 'limit' then we add the probeID for that bead type to a list which we later return.

Default is to use 14 as the limit as this is the 95th quantile for the chi-square distribution in question 

}
\value{
numeric vector giving probeIDs for all bead types on the array which give a chi-statistic > 'limit'

}

\author{Mark Dunning}


\examples{

data(BLData)

##Default setting is to find all bead types with a chi-statistic greater than 14

findHighestChis(BLData, 1)

##Now find all beads with chi-statistic for their distribution greater than 16

findHighestChis(BLData, 1, limit=16)

##ProbeIDs returned by the function can be examined in more detail

plotBeadLocations(BLData, probeID=278, array=1)

probeDiagnostics(BLData, 278, 1)


}

\keyword{methods}
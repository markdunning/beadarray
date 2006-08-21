\name{qsplineNormalise}
\alias{qsplineNormalise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Qsplinw Normalisation for BeadLevelList objects }
\description{
Normalizes expression intensities by fitting a smoothing curves to quantiles

Function is based on the 'normalize.qspline' function in 'affy'
}
\usage{
qsplineNormalise(BLData)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{can be either an BeadLevelList or BeadSummaryList object}
}
\details{
This normalization method uses the quantiles from each array and the target to fit a system of cubic splines to normalize the data. The target should be the mean (geometric) or median of each probe but could also be the name of a particular chip in the abatch object.

Parameters setting can be of much importance when using this method. The parameter fit.iter is used as a starting point to find a more appropriate value. Unfortunately the algorithm used do not converge in some cases. If this happens, the fit.iter value is used and a warning is thrown. Use of different settings for the parameter samples was reported to give good results. More specifically, for about 200 data points use samples = 0.33, for about 2000 data points use samples = 0.05, for about 10000 data points use samples = 0.02 (thanks to Paul Boutros).

The type argument should be one of "separate","pmonly","mmonly","together" which indicates whether to normalize only one probe type (PM,MM) or both together or separately.

}

\value{
Produces an BeadLevelList or BeadSummaryList (depending on input type) containing normalised intensity values. 

}
\author{Laurent and C Workman}


\examples{

data(BLData)

#BLData.q = qsplineNormalise(BLData)

#All distributions are fixed to be the same

#boxplot(log2(BLData.q$G)~col(BLData.q$G))

}

\references{Christopher Workman, Lars Juhl Jensen, Hanne Jarmer, Randy Berka, Laurent Gautier, Henrik Bj{o}rn Nielsen, Hans-Henrik Saxild, Claus Nielsen, S{o}ren Brunak, and Steen Knudsen. A new non-linear normal- ization method for reducing variability in dna microarray experiments. Genome Biology, accepted, 2002

}

\keyword{methods}


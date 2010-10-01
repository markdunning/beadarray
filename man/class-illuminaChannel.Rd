\name{illuminaChannel-class}
\docType{class}
\alias{illuminaChannel-class}

\title{Class "illuminaChannel" }

\description{
    A class to define how illumina bead-level data are summarized
}


\section{Slots/List Components}{
  Objects of this class contain the following slots
  \tabular{ll}{
  \code{transFun}:\tab A function that defines how to obtain a split list indexed by bead-type, containing all observations on a given array.\cr
  \code{outlierFun}:\tab A function that defines how.\cr
  \code{exprFun}:\tab A function for producing a single summary of expression level from a vector of bead-type intensities. e.g. mean .\cr
  \code{varFun}:\tab A function for producing a single summary of variability from a vector of bead-type intensities. e.g. sd \cr
  \code{name}: \tab Character vector that defines a name for the channel
  }
}

\author{Mark Dunning}

\examples{

greenChannelTransform = function(x){
	probeIDcol = which(colnames(x) == "ProbeID")
	gCol = which(colnames(x) == "G")

	split(x[,gCol], x[,probeIDcol])
}

redChannelTransform = function(x){
	probeIDcol = which(colnames(x) == "ProbeID")
	gCol = which(colnames(x) == "R")

	split(x[,gCol], x[,probeIDcol])
}

illuminaOutlierMethod = function(x){

M = median(x)

md = mad(x)

x[x < M+3*md & x > M-3*md]

}



greenChannel = new("illuminaChannel", greenChannelTransform, illuminaOutlierMethod, mean, sd, "G")
redChannel = new("illuminaChannel", redChannelTransform, illuminaOutlierMethod, mean, sd, "R")


###Example if we have two-channel data

logRatioTransform = function(x){
	probeIDcol = which(colnames(x) == "ProbeID")
	gCol = which(colnames(x) == "G")
	rCol = which(colnames(x) == "R")
	
	lr = x[,gCol]  - x[,rCol]
	
	split(x[,gCol],lr)
}


logRatio = new("illuminaChannel", logRatioTransform, illuminaOutlierMethod, mean, sd, "M")

}




\seealso{\code{\link{readIllumina}}}
\keyword{classes}

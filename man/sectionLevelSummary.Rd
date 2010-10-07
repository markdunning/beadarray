\name{sectionLevelSummary}
\alias{sectionLevelSummary}
\title{
    Create a summarized object
}
\description{
Function to summarize the data in a \code{beadLevelData} object into a form more ameanable for downstream analysis (with the same number of observations for each bead type).
}


\usage{
    sectionLevelSummary(BLData, channelList, probeIDs = NULL, weightNames ="wts", silent=FALSE)
}

\arguments{
    \item{BLData}{
        An object of class beadLevelData
    }
    \item{channelList}{
	List of objects of class illuminaChannel
   }

   \item{probeIDs}{
	Vector of ArrayAddressIDs to be included in the summarized object
	}

 	
   \item{weightNames}{
	name of column in the \code{beadLevelData} to take extract weights 
}

   \item{silent}{if FALSE the function will print information as it proceeds}	

}


\value{
    Returns an object of class ExpressionSetIllumina
}

\author{
    Mark Dunning
}

\examples{

	data(BLData)

	myMean = function(x) mean(x,na.rm=TRUE)
	mySd = function(x) sd(x,na.rm=TRUE)

	greenChannel = new("illuminaChannel", logGreenChannelTransform, illuminaOutlierMethod, myMean, mySd,"G")

	bsd = sectionLevelSummary(BLData, channelList = list(greenChannel))

	bsd


}
	


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{IO}


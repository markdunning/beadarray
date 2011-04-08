\name{checkPlatform}
\alias{checkPlatform}
\title{
Determine expression chip annotation
}
\description{

Works out which platform the bead level data represents given the probe IDs and compiled list of IDs from all expression platforms. This is based on the percentage of IDs on the array that overlap with IDs from known platforms. The platform with highest overlap is chosen. 

}


\usage{
    checkPlatform(BLData,verbose=FALSE)
}

\arguments{
    \item{BLData}{An object of class \code{\link{beadLevelData-class}}.
    \item{verbose}{If TRUE, report overlaps with known platforms}	
}
    
}

\value{
Character string of the best guess of what chip BLData is
}

\author{
Mark Dunning and Stewart MacArthur
}



\keyword{misc}


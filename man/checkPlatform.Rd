\name{checkPlatform}
\alias{checkPlatform}
\title{
Determine expression chip annotation
}
\description{

Works out which platform the bead level data represents given the probe IDs and compiled list of IDs from all expression platforms. 

}


\usage{
    checkPlatform(BLData)
}

\arguments{
    \item{BLData}{An object of class \code{\link{beadLevelData-class}}.

}
    
}

\value{
Character string of the best guess of what chip BLData is
}

\author{
Mark Dunning and Stewart MacArthur
}



\keyword{misc}


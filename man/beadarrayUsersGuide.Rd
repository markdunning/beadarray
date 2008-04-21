\name{beadarrayUsersGuide}
\alias{beadarrayUsersGuide}
\title{View beadarray User's Guide}
\description{
Finds the location of the beadarray User's Guide and opens it.
}
\usage{
beadarrayUsersGuide(view=TRUE, topic="beadlevel")
}
\arguments{
  \item{view}{logical, should the document be opened using the default PDF document reader? (default is \code{TRUE})}
  \item{topic}{character string specifying topic (\code{"beadlevel"} or \code{"beadsummary"})} 
}
\details{
  The function \code{vignette("beadarray")} will find the short beadarray Vignette which describes how to obtain the more detailed user's guide on the analysis of raw \code{"beadlevel"} or summarised \code{"beadsummary"} Illumina data.
   }
\value{Character string giving the file location.}

\seealso{
\code{\link[limma]{limmaUsersGuide}}
}

\author{Matt Ritchie}

\examples{
beadarrayUsersGuide(view=FALSE)
beadarrayUsersGuide(view=FALSE, topic="beadsummary")
}
\keyword{documentation}
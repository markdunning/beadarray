\name{backgroundCorrectSingleSection}
\alias{backgroundCorrectSingleSection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Background correct an array-section 
}
\description{
Function to perform background correction using a defined set of foreground and background intensities.
}
\usage{
backgroundCorrectSingleSection(BLData, array = 1, fg="Grn", bg="GrnB", newName = "Grn.bc")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{
a \code{beadLevelData} object
}
  \item{array}{
the number of the array-section to be corrected
}
  \item{fg}{
the name under which the foreground intensities are stored
}
  \item{bg}{
the name under which the background intensities are stored
}
  \item{newName}{
Name to store the corrected intensities 
}
}
\details{

After the foreground and background intensities are retrieved, the function calculates foreground - background

}
\value{
\code{beadLevelData} object with modified beadData slot for the particular section
}
\author{
Mark Dunning
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

##This will take the "Grn" and "GrnB" data for each section and save the corrected intensities in "Grn.bc"

data(BLData)

head(BLData[[1]])

for(i in 1:10){
	 BLData = backgroundCorrectSingleSection(BLData, array=i)
}


head(BLData[[1]])

}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

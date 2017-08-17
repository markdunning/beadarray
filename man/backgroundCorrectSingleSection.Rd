\name{backgroundCorrectSingleSection}
\alias{backgroundCorrectSingleSection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Background correct an array-section 
}
\description{
Function to perform a basic bead-level background correction using a defined set of foreground and background intensities.
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

This function takes two attributes of a bead-level object and returns that bead-level object with an additional attribute consisting of the difference of the other two. We anticipate this being used as a simple background correction step, returning the difference between foreground and background intensities to be used as the specific intensity associated with a bead. 

Functions to perform more advanced background correction steps can easily be constructed after consideration of this function.

}
\value{
\code{beadLevelData} object with modified beadData slot for the particular section
}
\author{
Mark Dunning
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{


\dontrun{

if(require(beadarrayExampleData)){

  data(exampleBLData)

  head(exampleBLData[[1]])

  for(i in 1:2){
	  exampleBLData = backgroundCorrectSingleSection(exampleBLData, array=i)
  }


  head(exampleBLData[[1]])

  } else {
  
  stop("You will need the beadarrayExampleData package to run this example")
}

}

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

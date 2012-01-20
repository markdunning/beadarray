\name{outlierplot}
\alias{outlierplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot outlier locations}
\description{
Function to plot where the outliers are located on a given array

}
\usage{
outlierplot(BLData, array = 1, transFun = logGreenChannelTransform, outlierFun = illuminaOutlierMethod,n=3, wtsname=NULL, horizontal = TRUE, nSegments = NULL, lowOutlierCol = "blue", highOutlierCol = "pink", outlierPch = ".", main = "", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{
a \code{beadLevelData} object
}
  \item{array}{
the number of the array to plot
}
  \item{transFun}{
a function defining how to transform the data prior to calculating outliers
}
  \item{outlierFun}{
function that will identify outliers
}
  \item{n}{
an indicator of how extreme an observation must be (e.g. how many MADs from the median), to be passed to the function that will identify outliers}
  \item{wtsname}{
column name of BLData object containing weights to feed to the outlier function
}
  \item{horizontal}{
if TRUE the longest edge of the array section will be on the x axis
}
  \item{nSegments}{
How many segments the section is divided into. If this argument is left as the default value (NULL) the code will attempt to extract this information from the relevant .sdf file.  If it can't be found then the segments will not be indicated on the final plot.
}
  \item{lowOutlierCol}{
what colour to plot outliers below the median
}
  \item{highOutlierCol}{
what colour to plot outliers above the median
}
  \item{outlierPch}{
plotting character for the outliers
}
  \item{main}{
an optional title for the plot
}
  \item{...}{
additional arguments
}


}
\details{

The function calls the specified outlier function to determine the outliers on the array and then plots their location. Points are coloured according the intensity of the bead is above or below the median for that bead-type. 

}
\value{
plot produced on current graphical device
}
\author{
Mark Dunning and Mike Smith
}

\examples{

if(require(beadarrayExampleData)){

data(exampleBLData)
outlierplot(exampleBLData, array=1, horizontal = FALSE)


}

}


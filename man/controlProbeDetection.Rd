\name{controlProbeDetection}
\alias{controlProbeDetection}

\title{
Percentage of beads detected
}
\description{
Function to calculate the percentage of beads matching a defined set of control types that are detected as having intensity above background level on an array-section.
}


\usage{
controlProbeDetection(BLData, transFun = logGreenChannelTransform, array = 1, controlProfile = NULL, tagsToDetect = list(housekeeping = "housekeeping", Biotin = "phage_lambda_genome", Hybridisation = "phage_lambda_genome:high"), negativeTag = "permuted_negative", detThresh = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{
a \code{beadLevelData} object
}
  \item{transFun}{
transformation to be applied to data
}
  \item{array}{
a numeric index of the array section
}
  \item{controlProfile}{
optional data frame defining ArrayAddressIDs belonging to each control type
}
  \item{tagsToDetect}{
vector of character strings defined which control types to interrogate
}
  \item{negativeTag}{
character string defining which control type to use as background
}
  \item{detThresh}{
numeric value for threshold for detection
}
}
\details{
Details of the controls on the array-section can be inferred from the annotation of the \code{beadLevelData} object or supplied as a data frame. The first column of the data frame should contain ArrayAddressIDs, with the control type of the each ID in the second column. The strings supplied in the \code{tagsToDetect} and \code{negativeTag} parameters should be present in this column. 

The ArrayAddressIDs that correspond to the specified tags are matching to the ArrayAddressIDs for the chosen array and intensities for all beads are extracted. The function implements Illumina's method for calculating the detection scores for all bead types on a given array. Within an array, Illumina discard negative control bead-types whose summary values are more than three MADs from the median for the negative controls. Illumina then rank the summarized intensity for each other bead-type against the summarized values for the remaining negative control bead-types and calculate a detection p-value  1-R/N, where R is the relative rank of the bead intensity when compared to the $N$ remaining negative controls. Thus, if a particular bead has higher intensity than all the negative controls it will be assigned a value of 0. This calculation is repeated for all arrays. 


The percentage reported is the percentage of beads of each control type that are detected at the defined threshold.

}

\author{
Mark Dunning
}

\seealso{\code{\link{beadStatusVector}},\code{\link{calculateDetection}}}


\examples{

if(require(beadarrayExampleData)){

data(exampleBLData)

for(i in 1:2){

print(controlProbeDetection(exampleBLData, array = i, tagsToDetect=c("housekeeping", "biotin"), negativeTag="negative"))

}

}
  
}


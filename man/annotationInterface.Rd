\name{getAnnotation}
\alias{annotation,beadLevelData-method}
\alias{annotation,ExpressionSetIllumina-method}
\alias{annotation<-,beadLevelData,character-method}
\alias{annotation<-,ExpressionSetIllumina,character-method}

\alias{getControlProfile}
\alias{suggestAnnotation}

\title{
Storage of annotation information
}

\usage{
    suggestAnnotation(data,verbose=FALSE)
}

\arguments{
    \item{data}{An object of class \code{\link{beadLevelData-class}}}.
    \item{verbose}{If TRUE, report overlaps with known platforms}	

} 


\description{
An interface to set or retrieve information about the annotation of a \code{beadLevelData} or \code{ExpressionSetIllumina} object. 
}

\details{

A character string is used to specify the annotation with the currently supported values being; Humanv4, Humanv3, Humanv2, Humanv1, Mousev2, Mousev1, Mousev1p1 and Ratv1. This string is used within beadarray to retrieve control probe IDs within particular QC functions.

The \code{suggestAnnotation} function tries to determine a suitable value of the annotation based on the probe IDs and compiled list of IDs from all expression platforms. This is based on the percentage of IDs on the array that overlap with IDs from known platforms. The platform with highest overlap is chosen. 

}

\examples{

if(require(beadarrayExampleData)){

data(exampleBLData)

annotation(exampleBLData)

suggestAnnotation(exampleBLData,verbose=TRUE)

annotation(exampleBLData) <- "Humanv2"
 

data(exampleSummaryData)

annotation(exampleSummaryData)

annotation(exampleBLData) <- "Humanv1"

}


}

\author{
Mark Dunning
}



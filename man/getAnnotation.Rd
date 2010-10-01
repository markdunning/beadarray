\name{getAnnotation}
\alias{getAnnotation}
\alias{setAnnotation}
\title{
Storage of annotation information
}
\description{
An interface to set or retrieve information about the annotation of a \code{beadLevelData} object. 
}
\usage{
getAnnotation(BLData)
setAnnotation(BLData, annoName)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{
a \code{beadLevelData} object
}
\item{annoName}{character string to set annotation}

}
\details{
A character string is used to specify the annotation with the currently supported values being; Humanv4, Humanv3, Humanv2, Humanv1, Mousev2, Mousev1, Mousev1p1 and Ratv1. This string is used within beadarray to retrieve control probe ids for each platform.
}
\value{
setAnnotation returns a \code{beadLevelData} object with the annotation stored.

getAnnotation returns the name of the annotation that is currently being stored.

}
\author{
Mark Dunning
}



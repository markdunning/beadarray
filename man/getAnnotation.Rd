\name{getAnnotation}
\alias{getAnnotation}
\alias{getControlAnno}
\alias{setAnnotation}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Storage of annotation information for Illumina expression chips}
\description{

Illumina use several control types for QA purposes, however the IDs of these controls change between different organisms and annotation revisions. Therefore we need to store the annotation of a chip in order to perform QA on the bead-level data.
The functions setAnnotation and getAnnotation are used to manage this annotation information.

}
\usage{
getAnnotation(BLData)
setAnnotation(BLData, aName)
}
\arguments{
  \item{BLData}{BeadLevelList for an Illumina expression array.}
  \item{aName}{Character to define the annotation of the chip}
}
\details{
 We currently store the annotation as a slot in a BeadLevelList. The value in the slot should match one of entries in \code{ExpressionControlData} (see example).}
\value{
setAnnotation returns a modified BeadLevelList with the new value for the annotation slot.
}
\author{Mark Dunning}

\examples{
data(BLData)
data(ExpressionControlData)
getAnnotation(BLData)
names(ExpressionControlData)

setAnnotation(BLData, "foo")
}

\seealso{\code{\link{ExpressionControlData}}}

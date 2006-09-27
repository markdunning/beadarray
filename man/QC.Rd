\name{QC}
\alias{QC}
\docType{data}
\title{Illumina quality control data}
\description{
An object of type \code{AssayData} containing quality control information for the example whole genome experiment.

The object has slots, Signal, StDev and Detection which are the average intensity, bead standard error and detection score for a variety of different controls placed on each array.
The column names of the Signal, StDev and Detection matrices are the names of each control used. See Illumina documentation for description of the controls and their interpretation.
}
\usage{data(BSData)}
\seealso{\code{\link[beadarray:readQC]{readQC}}}

\keyword{datasets}

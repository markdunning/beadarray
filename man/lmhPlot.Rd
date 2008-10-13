\name{lmhPlot}
\alias{lmhPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot the bead-level hybridisation controls}
\description{
Function for retrieving and plotting the hybridisation controls for an expression array. We know these controls should show high signal and are therefore useful for QA purposes. Moreover, we should expect to see a gradient between the low, medium and high controls. By considering all bead observations (unlike the plots produced by BeadStudio) we get an detailed impression of array quality. 
}
\usage{
lmhPlot(BLData, array = 1, plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{BeadLevelList object for an Illumina expression array which must have the annotation slot set appropriately.}
  \item{array}{The number of the array of interest}
  \item{plot}{If TRUE then a diagnostic plot will be produced, other only summary values will be returned.}
}
\details{
As well as plotting the controls, the function returns five measures that can be used to evaluate the quality of the array (see below). On a good quality array, we would expect to see 100\% for all these measures and a drop from 100\% could indicate a defect on the array. However, it should be noted that the HvsM percentage could drop below 100\% often due to the saturation effect often observed at high intensity. \\

Retrieving the control intensities relies on the accurate setting of the annotation slot of BLData, as the IDs of the controls change for different assays. If the wrong annotation is stored, we may not be able to find the controls, or misleading results may be produced. The function will try and give warning and error messages to guard against these problems. 
}

\value{
\item{LowDet - \%age of ``low'' control beads that are detected compared to the negative controls.}
\item{MedDet - \%age of ``medium'' control beads that are detected compared to the negative controls.}
\item{HighDet- \%age of ``high'' control beads that are detected compared to the negative controls.}
\item{MvsL - \%age of ``medium'' control beads that are detected compared to the ``low'' controls.}
\item{HvsM- \%age of ``high'' control beads that are detected compared to the ``medium'' controls.}
}

\references{see www.illumina.com/downloads/GX\_QualityControl\_TechNote.pdf for description of the hybridisation controls}
\author{Mark Dunning}
\seealso{\code{\link{calculateBeadLevelScores}},\code{\link{setAnnotation}} }

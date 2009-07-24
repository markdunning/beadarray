\name{calculateBeadLevelScores}
\alias{calculateBeadLevelScores}
\alias{gradientPlot}

\title{Quality assessment for expression chips }
\description{
A collection of functions for tabulating and plotting various quality control measurements derived from the bead-level data for Illumina expression chips. Currently, Humanv1, Humanv2, Humanv3, Mousev1, Mousev1p1, Mousev2 and Rat chips are supported.
}
\usage{
calculateBeadLevelScores(BLData, path = "QC", log = TRUE, plot = FALSE, 
                         replacePlots = TRUE, writeToFile = TRUE,
                         fileType = c(".jpeg", ".pdf"))
}
\arguments{
  \item{BLData}{A BeadLeveList for an expression chip. The annotation slot of the object should define the type of chip }
  \item{path}{Specifies the directory where diagnostics plot are to be saved in}
  \item{log}{(used for outlier calculations) if TRUE calculate outliers on the log2 scale.  If FALSE calculate outliers on the original scale}
  \item{plot}{if TRUE then diagnostic plots will be generated}
  \item{writeToFile}{Argument describing whether results of QA assessment should be output to html, txt or not output at all}
  \item{replacePlots}{if TRUE any plots that have already been created will be replaced}
  \item{fileType}{File format for images to be written to. Currently either pdf or jpeg}
}
\details{For these QA tools we make use of the controls probes that Illumina
use on their expression chips to detect the presence, or lack of,
expression. See \url{www.illumina.com/downloads/GX_QualityControl_TechNote.pdf} for an overview of these controls. Illumina provide a means to visualize these controls, but the values reported are \emph{after outlier removal} and there is no way to infer how many outliers are removed. Therefore, one does not get a true impression of the quality of an array. For instance, low intensity observations for positive controls may indicate a spatial defect.

For our QA measurements we perform a detection score calculation the same as Illumina, except on per-bead observations for each control type rather than the summarized values. Specifically, we test each bead observation of a given control bead-type for detection by computing a p-value: \eqn{1-R/N}, where \eqn{R} is the relative rank of the bead intensity when compared to the \eqn{N} negative controls. Thus, if a particular bead has higher intensity than all the negative controls it will be assigned a value of \eqn{0}. After these p-values have been calculated for all replicates of the bead type we report the percentage of beads with p-values lower than a set threshold of 0.05 (currently in favour in the Illumina literature). The percentage of beads that are detected at a set threshold is then reported. Another adaptation is to change the bead-types used as a reference in the calculation rather than the negative controls. For example, there are a series of sample-independent controls that have probe sequences complementary to oligonucleotides spiked into the hybridization solution and hence should always have detectable signal. For some of these bead types (six on the Human6 V3), the concentration is either "medium", "low" or "high", with the intention that there should be a predictable gradient between the controls. Thus, we test if the bead-types with a medium concentration are detected compared to the low controls and similar for the medium and high controls.

The purpose of calculateBeadLevelScores is to calculate the following QA measures for all arrays in the BeadLevelList objects and return them in the arrayInfo slot of the BeadLevelList object. We also record the number of outliers found on the array.

If the plot argument to calculateBeadLevelScores is set to TRUE, then a number of diagnostic plots will be produced for each and compiled into a HTML page for that array. The location of these completed pages is specified by the path argument. Finally, if writeToFile is set to html, a html page compiling all the QA measures for the chip will be created.

}
\value{
A modified version of BeadLevelList is created with the QA measures stored in
a qcScores slot.

\item{HkpDet}{\%age of housekeeping control beads that are detected compared to the negative controls.}
\item{BioDet}{\%age of  biotin labelling control beads that are detected compared to the negative controls.}
\item{LowDet}{\%age of "low" control beads that are detected compared to the negative controls.}
\item{MedDet}{\%age of "medium" control beads that are detected compared to the negative controls.}
\item{HighDet}{\%age of "high" control beads that are detected compared to the negative controls.}
\item{MvsL}{\%age of "medium" control beads that are detected compared to the "low" controls.}
\item{HvsM}{\%age of "high" control beads that are detected compared to the "medium" controls.}
}

\author{Mark Dunning and Andy Lynch}
\seealso{\code{\link{outlierPlot}}, \code{\link{lmhPlot}}, \code{\link{poscontPlot}}, \code{\link{backgroundControlPlot}}}

\name{poscontPlot}
\alias{poscontPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot bead-level housekeeping and biotin controls}
\description{
Function for retrieving and plotting the biotin and housekeeping controls for an expression array. We know these controls should show high signal and are therefore useful for QA purposes. The housekeeping control targets a bead-type believed to be universally expressed whereas the biotin control targets the biotin used for staining. By considering all bead observations (unlike the plots produced by BeadStudio) we get an detailed impression of array quality. 
}
\usage{
poscontPlot(BLData, array = 1, plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{A BeadLevelList object for an Illumina expression array}
  \item{array}{The number of the array of interest}
  \item{plot}{if TRUE a diagnostic plot will be produced, otherwise only summary values are returned.}
}
\details{

The IDs for the housekeeping and biotin controls are retrieved making use of the annotation information stored for the array. If this information is incorrect, or missing, then the function will not proceed correctly. Valid names for the annotation slot must match the names of the ExpressionControlData object provided with beadarray. Assuming the controls can be identified, by then perform a detection test for each bead observation by computing a p-value:- $1-R/N$, where $R$ is the relative rank of the bead intensity when compared to the $N$ negative controls. Thus, if a particular bead has higher intensity than all the negative controls it will be assigned a value of $0$. After these p-values have been calculated for all replicates of the bead type we report the percentage of beads with p-values lower than a set threshold of 0.05 (currently in favour in the Illumina literature). The percentage of beads that are detected at a set threshold is then reported for the biotin housekeeping and biotin controls respecitvely.\\

If plot is TRUE, the values of all bead observations are plotted for each bead-type respectively. Any beads with low intensity may be due to array defects and would be impossible to detect with summarized data only.\\

An important point to note that the utility of these housekeeping controls is dependant on the probe sequences used targeting the intended genes. If a particular housekeeping control shows lower expression, then it is worth checking the probe sequence by a BLAT search of by checking pre-complied tables from {\tt http://www.compbio.group.cam.ac.uk/Resources/Annotation/index.html}.

}

\value{
\begin{itemize}
\item{HkpDet}{ - \%age of housekeeping control beads that are detected compared to the negative controls.}
\item{BioDet} {- \%age of  biotin labelling control beads that are detected compared to the negative controls.}
\end{itemize}
}

\references{ ~put references to the literature/web site here ~ }
\author{Mark Dunning and Andy Lynch}
\seealso{ \code{\link{calculateBeadLevelScores}}, \code{\link{setAnnotation}}}


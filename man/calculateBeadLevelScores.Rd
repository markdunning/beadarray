\name{calculateBeadLevelScores}
\alias{calculateBeadLevelScores}
\alias{probePairsPlot}
\alias{lmhPlot}
\alias{backgroundControlPlot}
\alias{poscontPlot}
\alias{gradientplot}
\alias{outlierplot}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Quality assessment for expression chips }
\description{
A collection of functions for tabulating and plotting various quality control measurements derived from the bead-level data for Illumina expression chips. Currently, HumanV2, HumanV3, MouseV1 and MouseV2 chips are supported
}
\usage{
calculateBeadLevelScores(BLData, path = "QC", plot = FALSE, writeToFile = c(NULL, "html", "txt"))
probePairsPlot(BLData,array=1,plot=FALSE)
lmhPlot(BLData,array=1,plot=FALSE)
backgroundControlPlot(BLData,array=1,plot=FALSE)
poscontPlot(BLData,array=1,plot=FALSE)
gradientplot(BLData,array=array,plot=FALSE)
outlierplot(BLData,array=array,plot=FALSE)

}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{A BeadLeveList for an expression chip. The annotation slot of the object should define the type of chip }
  \item{path}{Specifies the directory where diagnostics plot are to be saved in}
  \item{plot}{if TRUE then diagnostic plots will be generated}
  \item{writeToFile}{ ~~Describe \code{writeToFile} here~~ }
}
\details{The QA measures calculated by calculateBeadLevelScores are:-

\begin{itemize}
\item[HkpDet]{ - \%age of housekeeping control beads that are detected compared to the negative controls.}
\item[BioDet]{ - \%age of  biotin labelling control beads that are detected compared to the negative controls.}
\item[LowDet]{ - \%age of ``low'' control beads that are detected compared to the negative controls.}
\item[MedDet]{ - \%age of ``medium'' control beads that are detected compared to the negative controls.}
\item[HighDet]{ - \%age of ``high'' control beads that are detected compared to the negative controls.}
\item[MvsL]{ - \%age of ``medium'' control beads that are detected compared to the ``low'' controls.}
\item[HvsM]{ - \%age of ``high'' control beads that are detected compared to the ``medium'' controls.}
\end{itemize}

These are stored in the qcScores slot of the BeadLevelList object

}
\value{
A modified version of BeadLevelList is created
}
\references{ ~put references to the literature/web site here ~ }
\author{Mark Dunning and Andy Lynch}


\name{readQC}
\alias{readQC}


\title{Read Illumina control intensities}
\description{
Reads the standard format of Illumina control intensities output by BeadStudio
}

\usage{
readQC(file, sep="\t", skip=8, controlID = "ProbeID", columns = list(exprs = "AVG_Signal", se.exprs="BEAD_STDERR", NoBeads = "Avg_NBEADS", Detection="Detection Pval"), dec=".", quote="")

}
\arguments{
  \item{file}{character string giving the name of the file output by BeadStudio
    containing the control probe intensities.  This file should be either the 'ControlProbeProfile' or 'ControlGeneProfile'.} 
  \item{sep}{a character string for the file separator}
  \item{skip}{number of lines of header information to ignore in the file}
  \item{controlID}{character string specifying the column that contains the 
    (unique) control probe IDs}
  \item{columns}{a vector of column names to read from the file} 
  \item{dec}{the character used in the file for decimal points}  
  \item{quote}{the set of quoting characters (disabled by default)}
}
\details{

  The format of the quality control files differs slightly between BeadStudio 
  versions 1 and later versions.  This function is able to read in data 
  in either format

  Note that if the control identifiers are non-unique, the duplicated 
  rows are removed.  This may occur if the 'ControlProbeProfile' is 
  exported from BeadStudio and \code{controlID="TargetID"} is specified 
  (the "ProbeID" column has a unique identifier in the 'ControlProbeProfile',
  whereas the "TargetID" may not, as multiple beads can be of the same type).

  Once read in, the control intensities can be used for quality assessment purposes.
}

\value{
  \code{readQC} produces an \code{assayData} object with a list of items defined 
  by the \code{columns} parameter.  The row names of each matrix are given by the 
  \code{controlID} argument .
}
\author{Mark Dunning}

\examples{
##Code to read the example quality control file included with the
#package.
#QC = readQC("ControlGeneProfile.txt", controlID="TargetID")
#the average expression of each control can then be accessed by the $ operator
#QC$exprs
}

\keyword{methods}

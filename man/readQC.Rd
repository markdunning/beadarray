\name{readQC}
\alias{readQC}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read Illumina quality control information}
\description{
Reads the standard format of Illumina quality control information and produces diagnostic plots
}

\usage{
readQC(file, columns = list(exprs = "AVG_Signal", NoBeads = "Avg_NBEADS", Detection="Detection", BeadStDev="BEAD_STDERR", Narrays="NARRAYS", arrayStDev = "ARRAY_STDEV", controlID = "ProbeID", controlType="TargetID"),sep="\t",skip=7,header=T)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{name of file containing qc information} 
  \item{columns}{a vector of column names to read from the file} 
  \item{skip}{number of lines of header information to ignore in the file}
  \item{header}{if TRUE the column names in the file are read}
  \item{sep}{a character string for the file separator}
		 
}
\details{

The format of the quality control
  files differs slightly between BeadStudio versions 1 and 2. Version 1 of the
  software gives one averaged value for each control type, whereas version 2
  gives suumarised values for each control of a particular type. The user does
  not have to know the version of BeadStudio used to generate the file.


  Once read in, the quality control information can be very useful for
  diagnostic purposes.
}

\value{
  readQC produces an assayData object with list items defined by the
  columns parameter.

 
}
\author{Mark Dunning}

\examples{
##Code to read the example quality control file included with the
#package.
#QC = readQC(sep=",", skip=7, columns=list(exprs="AVG.Signal", BeadStDev="SeqVAR"))

#the average expression of each control can then be accessed by the $
operator

#QC$exprs


}

\keyword{methods}
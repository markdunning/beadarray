\name{readQC}
\alias{readQC}


\title{Read Illumina quality control information}
\description{
Reads the standard format of Illumina quality control information and produces diagnostic plots
}

\usage{
readQC(file, sep="\t", skip=8, controlID = "ProbeID", columns = list(exprs = "AVG_Signal", se.exprs="BEAD_STDERR", NoBeads = "Avg_NBEADS", Detection="Detection Pval"), dec=".", quote="")

}
\arguments{
  \item{file}{name of file containing qc information} 
  \item{sep}{a character string for the file separator}
  \item{skip}{number of lines of header information to ignore in the file}
  \item{controlID}{character string specifying the column that contains the 
    (unique) control probe IDs}
  \item{columns}{a vector of column names to read from the file} 
  \item{dec}{the character used in the file for decimal points}  
  \item{quote}{the set of quoting characters (disabled by default)}
}
\details{

  The format of the quality control
  files differs slightly between BeadStudio versions 1 and 2. Version 1 of the
  software gives one averaged value for each control type, whereas version 2
  gives summarised values for each control of a particular type. The user does
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
#QC = readQC("raw_data_qcinfo.csv", sep=",", skip=7, columns=list(exprs="AVG.Signal", se.exprs="SeqVAR"))

#the average expression of each control can then be accessed by the $ operator

#QC$exprs


}

\keyword{methods}

\name{readBeadSummaryData}
\alias{readBeadSummaryData}
\title{Read BeadStudio gene expression output}
\description{
Function to read the output of Illumina's BeadStudio software into beadarray
}
\usage{
readBeadSummaryData(dataFile, qcFile=NULL, sampleSheet=NULL,
                    sep="\t", skip=8, ProbeID="ProbeID",
                    columns = list(exprs = "AVG_Signal", se.exprs="BEAD_STDERR",
                        NoBeads = "Avg_NBEADS", Detection="Detection Pval"),
                    qc.sep="\t", qc.skip=8, controlID="ProbeID", 
                    qc.columns = list(exprs="AVG_Signal", se.exprs="BEAD_STDERR", 
			NoBeads="Avg_NBEADS", Detection="Detection Pval"), 
		    annoPkg=NULL, dec=".", quote="")

}
\arguments{
  \item{dataFile}{character string specifying the name of the file containing
    BeadStudio output}
  \item{qcFile}{character string giving the name of the file containing the 
    control probe intensities (optional)}
  \item{sampleSheet}{character string used to specify the file containing sample
    infomation (optional)}
  \item{sep}{field separator character for the dataFile (\code{"\t"} for tab delimited or \code{","} for comma separated)}
  \item{skip}{number of header lines to ignore at the top of \code{dataFile}}
  \item{ProbeID}{character string of the column in \code{dataFile} that contains identifiers 
    that can be used to uniquely identify each probe}
  \item{columns}{list defining the column headings in \code{dataFile} which 
    correspond to the matrices stored in the \code{assayData} slot of the final \code{ExpressionSetIllumina} object}
  \item{qc.sep}{separator character for qcFile}
  \item{qc.skip}{number of lines to skip at the top of qcFile}
  \item{controlID}{character string specifying the column in \code{qcFile} that contains 
     the identifiers that uniquely identify each control probe}
  \item{qc.columns}{list defining the column headings in \code{qcFile} which 
     correspond to the matrices stored in the \code{QCInfo} slot of
     the final \code{ExpressionSetIllumina} object}
  \item{annoPkg}{character string specifying the name of the annotation package 
     (only available for certain expression arrays at present)}
  \item{dec}{the character used in the file for decimal points}
  \item{quote}{the set of quoting characters (disabled by default)}
  }

\details{
  This function can be used to read gene expression data exported
  from versions 1,2 and 3 of the Illumina BeadStudio application.
  The format of the BeadStudio output may depend on the version number.
  For example, the file may be comma or tab seperated of have header
  information at the top of the file. The parameters sep and skip
  can be used to adapt the function as required (i.e. skip=7 is 
  appropriate for data from earlier version of BeadStudio, and skip=0 is
  required if header information hasn't been exported.

  The format of the BeadStudio file is assumed to have one row for each
  probe sequence in the experiment and a set number of columns for each
  array. The columns for each array are chosen from the user within
  BeadStudio. These can include the average expression for each probe,
  the standard error and the number of beads and more.

  It is assumed that the average bead intensities for each
  array appear in columns with headings of the form 'AVG\_Signal-ARRAY1',
  'AVG\_Signal-ARRAY2',...,'AVG\_Signal-ARRAYN' for the N arrays found in the
  file. All other columns headings are found in the same way.

  NOTE:  With version 2 of BeadStudio it is possible to export annotation and
  sequence information along with the probe intensities.  We \_don't\_ recommend 
  exporting this information, as special characters found in the annotation 
  columns can cause problems when reading in the data.  This annotation information
  can be retrieved later on from other Bioconductor packages.

  The default object created by readBeadSummaryData is an
  ExpressionSetIllumina object. The columns parameter can be used to select 
  the relevant column headings to read in.
  
  If quality control information has been exported from BeadStudio
  (control probe profile) this may be read into beadarray too. The
  qc.skip and qc.sep and qc.columns parameters can be used to adjust
  for the contents of the file.  If the 'ControlGeneProfile' is exported, 
  from BeadStudio, you will need to set controlID="TargetID" in qc.columns.

  Sample sheet information can also be used. This is a file format used
  by Illumina to specify the contents of each array in the experiment.

  Note that if the probe identifiers are non-unique, the duplicated 
  rows are removed.  This may occur if the 'SampleGeneProfile' is 
  exported from BeadStudio and/or ProbeID="TargetID" is specified 
  (the "ProbeID" column has a unique identifier in the SampleProbeProfile,
  whereas the "TargetID" may not, as multiple beads can target the same 
  transcript).
}

\value{
  An ExpressionSetIllumina object with the slots defined in the columns
  parameter filled with data from the input file.

}

\author{Mark Dunning and Mike Smith}

\seealso{
  \code{\link[beadarray:ExpressionSetIllumina]{ExpressionSetIllumina}}
}
\examples{
##code to read the example BeadStudio (version 2) output distributed with the package
#dataFile = "raw_data.csv"
#sampleSheet = "raw_data_sample_sheet.csv"
#qcFile = "raw_data_qcinfo.csv"
#BSData <- readBeadSummaryData(dataFile, qcFile=qcFile, sampleSheet=sampleSheet, skip=7, columns=list(exprs="AVG_Signal", se.exprs="BEAD_STDEV", NoBeads="Avg_NBEADS"), qc.columns=list(exprs="AVG.Signal", se.exprs="SeqVAR"),qc.sep="," ,sep=",", qc.skip=7)
}
\keyword{IO}

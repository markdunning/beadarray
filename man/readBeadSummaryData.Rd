\name{readBeadSummaryData}
\alias{readBeadSummaryData}
\title{Read BeadStudio gene expression output}
\description{
Function to read the output of Illumina's software into beadarray
}
\usage{
readBeadSummaryData(dataFile, qcFile=NULL, sampleSheet=NULL,
header=TRUE, sep="\t", ProbeID="TargetID",skip=8,
                                   columns = list(exprs = "AVG_Signal", NoBeads = "Avg_NBEADS",
                                     Detection="Detection", se.exprs="BEAD_STDERR", Narrays="NARRAYS", arrayStDev = "ARRAY_STDEV"),
                                   qc.columns = list(ControlID="TargetID",exprs="AVG_Signal", Detection="Detection", Narrays="NARRAYS", se.exprs="BEAD_STDERR", NoBeads="Avg_NBEADS", arrayStDev="ARRAY_STDEV"), annoPkg=NULL,
                                   qc.sep="\t", qc.skip=8, ...)

}
\arguments{
  \item{dataFile}{character string specifying the name of the file containing
    BeadStudio output}
  \item{qcFile}{optional character string giving the location of the
    control probe profile exported from BeadStudio}
  \item{sampleSheet}{optional character string used to specify a sample
    infomation file}
  \item{sep}{seperator character for the dataFile}
  \item{header}{set to TRUE if the dataFile contains column headings}
  \item{skip}{number of header lines to ignore at the top of dataFile}
  \item{ProbeID}{name of the column heading that contains the probe IDs
    that can be used to uniquely identify each row}
  \item{columns}{list defining which column headings are to be read from
    dataFile and which slots in the ExpressionSetIllumina object they
    correspond to}
  \item{qc.columns}{list defining which column headings are to be read
    from the qcFile and which items in the QCInfo slot of
     ExpressionSetIllumina they correspond to}  
  \item{annoPkg}{character string specifying the annotation package 
     name for the arrays being read in (only available for certain 
     expression arrays at present)}
  \item{qc.sep}{seperator character for qcFile}
  \item{qc.skip}{number of lines to skip at the top of qcFile}
  \item{...}{further arguments to be passed to \code{read.table}}
  }

\details{ 
  This function can be used to read gene expression data exported
  from versions 1,2 and 3 of the Illumina BeadStudio application.
  The format of the BeadStudio output may depend on the version number.
  For example, the file may be comma or tab seperated of have header
  information at the top of the file. The parameters sep and skip
  can be used to adapt the function as required.  Currently, the default
  options have been set to handle version 3 output.

  The format of the BeadStudio file is assumed to have one row for each
  probe sequence in the experiment and a set number of columns for each
  array. The columns for each array are chosen from the user within
  BeadStudio. These can include the average expression for each probe,
  the standard error and the number of beads and more.

  With version 2 of BeadStudio it is possible to export annotation and
  sequence information. We recommend this data is not exported if the
  file is to be read into beadarray as it can be retrieved later on from
  other Bioconductor packages.

  It is assumed that the average bead intensities for each
  array appear in columns with headings of the form 'AVG\_Signal-ARRAY1',
  'AVG\_Signal-ARRAY2',...,'AVG\_Signal-ARRAYN' for the N arrays found in the
  file. All other columns headings are found in the same way.

  The default object created by readBeadSummaryData is an
  ExpressionSetIllumina object with enough slots to store all the
  columns that it is possible to export from BeadStudio. The columns
  parameter can be used to select which slots can be found from the data
  and the relevant column headings.
  
  If quality control information has been exported from BeadStudio
  (control probe profile) this may be read into beadarray too. The
  qc.skip and qc.sep and qc.columns parameters can be used to adjust
  for the contents of the file.
  
  Sample sheet information can also be used. This is a file format used
  by Illumina to specify the contents of each array in the experiment.
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
#BSData <- readBeadSummaryData(dataFile, qcFile=qcFile, sampleSheet=sampleSheet,skip=7, columns=list(exprs="AVG_Signal", se.exprs="BEAD_STDEV", NoBeads="Avg_NBEADS"), qc.columns=list(exprs="AVG.Signal", se.exprs="SeqVAR"),qc.sep="," ,sep=",", qc.skip=7)
}
\keyword{IO}
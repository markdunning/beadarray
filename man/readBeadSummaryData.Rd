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
                        nObservations = "Avg_NBEADS", Detection="Detection Pval"),
                    qc.sep="\t", qc.skip=8, controlID="ProbeID", 
                    qc.columns = list(exprs="AVG_Signal", se.exprs="BEAD_STDERR", 
			nObservations="Avg_NBEADS", Detection="Detection Pval"), 
		    illuminaAnnotation=NULL, dec=".", quote="", annoCols = c("TargetID", "PROBE_ID","SYMBOL"))

}
\arguments{
  \item{dataFile}{character string specifying the name of the file containing the 
    BeadStudio output for each probe on each array in an experiment (required).  
    Ideally this should be the 'SampleProbeProfile' from BeadStudio.}
  \item{qcFile}{character string giving the name of the file containing the 
    control probe intensities (optional).  This file should be either the 
    'ControlProbeProfile' or 'ControlGeneProfile' from BeadStudio.}
  \item{sampleSheet}{character string used to specify the file containing sample
    infomation (optional)}
  \item{sep}{field separator character for the \code{dataFile} (\code{"\t"} for 
   tab delimited or \code{","} for comma separated)}
  \item{skip}{number of header lines to skip at the top of \code{dataFile}.  
   Default value is 8.}
  \item{ProbeID}{character string of the column in \code{dataFile} that contains 
    identifiers that can be used to uniquely identify each probe}
  \item{columns}{list defining the column headings in \code{dataFile} which 
    correspond to the matrices stored in the \code{assayData} slot of the final \code{ExpressionSetIllumina} object}
  \item{qc.sep}{field separator character for \code{qcFile}}
  \item{qc.skip}{number of header lines to skip at the top of \code{qcFile}}
  \item{controlID}{character string specifying the column in \code{qcFile} that contains 
     the identifiers that uniquely identify each control probe}
  \item{qc.columns}{list defining the column headings in \code{qcFile} which 
     correspond to the matrices stored in the \code{QCInfo} slot of
     the final \code{ExpressionSetIllumina} object}
  \item{illuminaAnnotation}{character string specifying the name of the annotation package 
     (only available for certain expression arrays at present)}
  \item{dec}{the character used in the \code{dataFile} and \code{qcFile} for decimal points}
  \item{quote}{the set of quoting characters (disabled by default)}
  \item{annoCols}{additional columns containing annotation to be read from the file}


  }

\details{
  This function can be used to read gene expression data exported
  from versions 1,2 and 3 of the Illumina BeadStudio application.
  The format of the BeadStudio output will depend on the version number.
  For example, the file may be comma or tab separated of have header
  information at the top of the file. The parameters \code{sep} and \code{skip}
  can be used to adapt the function as required (i.e. skip=7 is 
  appropriate for data from earlier version of BeadStudio, and skip=0 is
  required if header information hasn't been exported.

  The format of the BeadStudio file is assumed to have one row for each
  probe sequence in the experiment and a set number of columns for each
  array. The columns which are exported for each array are chosen by the 
  user when running BeadStudio.  At a minimum, columns for average intensity
  standard error, the number of beads and detection scores should be exported, 
  along with a column which contains a unique identifier for each bead type 
  (usually named "ProbeID").

  It is assumed that the average bead intensities for each array appear in 
  columns with headings of the form 'AVG\_Signal-ARRAY1',
  'AVG\_Signal-ARRAY2',...,'AVG\_Signal-ARRAYN' for the N arrays found in the
  file.  All other column headings are matched in the same way using the character 
  strings specified in the \code{columns} argument.

  NOTE:  With version 2 of BeadStudio it is possible to export annotation and
  sequence information along with the intensities.  We \_don't\_ recommend 
  exporting this information, as special characters found in the annotation 
  columns can cause problems when reading in the data.  This annotation information
  can be retrieved later on from other Bioconductor packages.

  The default object created by readBeadSummaryData is an
  \code{ExpressionSetIllumina} object.
  
  If the control intensities have been exported from BeadStudio
  ('ControlProbeProfile') this may be read into beadarray as well. The
  \code{qc.skip}, \code{qc.sep} and \code{qc.columns} parameters can be 
  used to adjust for the contents of the file.  If the 'ControlGeneProfile' 
  is exported, you will need to set \code{controlID="TargetID"}.

  Sample sheet information can also be used. This is a file format used
  by Illumina to specify which sample has been hybridised to each array 
  in the experiment.

  Note that if the probe identifiers are non-unique, the duplicated 
  rows are removed.  This may occur if the 'SampleGeneProfile' is 
  exported from BeadStudio and/or \code{ProbeID="TargetID"} is specified 
  (the "ProbeID" column has a unique identifier in the 'SampleProbeProfile',
  whereas the "TargetID" may not, as multiple beads can target the same 
  transcript).
}

\value{
  An \code{ExpressionSetIllumina} object.
}

\author{Mark Dunning and Mike Smith}

\examples{
##Read the example data from
##http://www.switchtoi.com/datasets/asuragenmadqc/AsuragenMAQC_BeadStudioOutput.zip
##To follow this example, download the zip file 


\dontrun{
dataFile = "AsuragenMAQC-probe-raw.txt"

qcFile = "AsuragenMAQC-controls.txt"

BSData = readBeadSummaryData(dataFile=dataFile, qcFile=qcFile, controlID="ProbeID",skip=0,qc.skip=0, qc.columns=list(exprs = "AVG_Signal"))

}

}
\keyword{IO}

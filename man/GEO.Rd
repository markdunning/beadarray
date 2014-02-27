\name{makeGEOSubmissionFiles}
\alias{makeGEOSubmissionFiles}
\alias{createGEOMeta}
\alias{createGEOMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create files for a Gene Expression Omnibus submission
}
\description{
Following the guidelines at http://www.ncbi.nlm.nih.gov/geo/info/geo_illu.html, we provide functionality to create the three types of file required for a GEO submission of Illumina expression data. These are 1) sample metadata 2) Normalised intensities 3) Raw intensities. 


}
\usage{
makeGEOSubmissionFiles(normData, rawData, forceDetection = TRUE, basename = "GEO", softTemplate = NULL)
createGEOMeta(normData,basename)
createGEOMatrix(normData,forceDetection=TRUE,basename="GEO",softTemplate=NULL,normalised=T)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{normData}{
A normalised ExpressionSetIllumina object
}
  \item{rawData}{
A 'raw' non-normalised ExpressionSetIllumina object
}
  \item{forceDetection}{
TRUE forces the calculation of detection scores (if they are not present).
}
  \item{basename}{
A prefix that all the output files will be given
}
  \item{softTemplate}{
A SOFT file downloaded from GEO that we wish to annotate against. e.g. GPL10558.annot for Humanv4
}
\item{normalised}{
if TRUE the suffix ProcessedMatrix will be used to write the intensities to. Otherwise the Raw suffix is used.
}
}
\details{
Each file is written as a tab-delimited file in the current working directory. Note that the metadata contains the correct number of columns for the number of samples in the data, but the contents should be entered / curated manually prior to submission.\\

Meta file: A template file that can be used to record all the relevant metadata for the experiment, plus the information required by GEO. The SAMPLES section of the file contains the per-sample metadata and has columns generate from the required GEO data an the phenoData stored with the input data\\

Processed Matrix - GEO seems to want normalised intensities and detection scores for each probe on each array. Detection scores can be generated if they are not present in the data. The probes that are reported in this file are decided using the annotation slot of the object and querying the appropriate annotation package. Otherwise, a SOFT file downloaded from GEO can be used. \\
Raw Matrix - As above, but non-normalised intensities are used

}
\value{
}
\references{
http://www.ncbi.nlm.nih.gov/geo/info/geo_illu.html
}
\author{
Mark Dunning
}
\examples{
\dontrun{
if(require(beadarrayExampleData)){
data(exampleSummaryData)
rawdata <- channel(exampleSummaryData, "G")
normdata <- normaliseIllumina(rawdata)
makeGEOSubmissionFiles(normdata, rawdata)


##To annotate against a SOFT file that can be downloaded from GEO
##e.g. ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz

makeGEOSubmissionFiles(normdata, rawdata,softTemplate="GPL10558.annot")
}

}

}
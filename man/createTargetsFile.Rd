\name{createTargetsFile}
\alias{createTargetsFile}

\title{
A function to generate a targets file given a directory of Illumina bead-level files
}
\description{
This function, when pointed to a directory containing Illumina bead-level files (e.g. txt, idat, locs, tif) will return a simple targets file of the sort expected by beadarray. Note that a user created targets file is likely to be of greater value.
}
\usage{
createTargetsFile(dir = NULL, nochannels = 1, channel1 = "Grn", channel2 = "Red", txtsuff = "txt", imgsuff = "tif", locssuff = "locs", xmlsuff = "xml", verbose = FALSE, special = c("sdf", "fiducial"), ColourConfusionStop = T, metricsflag = "Metrics", metsep = "\t", metricsection = "Section", metricchip = "Matrix")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dir}{
\code{dir}: The directory containing the Illumina bead-level files. By default, will search the working directory.
}
  \item{nochannels}{
\code{nochannels}: Does the directory contain 1 or 2 channel arrays? Setting this argument to be null will result in the function making its best guess.
}
  \item{channel1}{
\code{channel1}: The string indicating that files are associated with the first channel (usually Grn). 
}
  \item{channel2}{
\code{channel2}: The string indicating that files are associated with the second channel (usually Red).  
}
  \item{txtsuff}{
\code{txtsuff}: The suffix of files containing the bead-level intensities (usually txt, but occasionally csv).
}
  \item{imgsuff}{
\code{imgsuff}: The suffix of files containing the images.
}
  \item{locssuff}{
\code{locssuff}: The suffix of files containing the precise bead locations (usually locs).
}
  \item{xmlsuff}{
\code{xmlsuff}: The suffix of files containing the meta-data (usually xml).
}
  \item{verbose}{
\code{verbose}: Determines whether or not the function reports on its progress as it goes along.
}
  \item{special}{
\code{special}: Files with names containing special words (such as fiducial) are ignored.
}
  \item{ColourConfusionStop}{
\code{ColourConfusionStop}: This determines the behaviour of the function if there is a discrepancy between the number of channels specified, and the number apparently present. 
}
  \item{metricsflag}{
code{metricsflag}: This gives the key word that can be used to identify metrics files.  
}
  \item{metsep}{
\code{metsep}: This gives the cell separator used in the metrics file.   
}
  \item{metricsection}{
\code{metricsection}: This gives the column heading used in the metrics file to indicate array section names.   
}
  \item{metricchip}{
\code{metricchip}: This gives the column heading used in the metrics file to indicate the chip name. 
}
}
\details{
This function bases its resultant targets file on the files with suffix txtsuff.
}
\value{
This returns a dataframe containing 
}

\author{
Andy Lynch
}


\seealso{
readIlluminaData()
}
\examples{

#createTargetsFile(verbose=T)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }


\name{readBGX}
\alias{readBGX}

\title{Read Illumina .bgx file into R}
\description{
Reads in an unzipped Illumina .bgx file, which provides further information
on each bead type, including the controls.
}
\usage{
readBGX(filename, path=".", sep="\t", quote="", header=TRUE,
           probeStart="\\\[Probes\\\]", controlStart="\\\[Controls\\\]", ...)
}
\arguments{
  \item{filename}{character vector specifying name of unzipped bgx file}
  \item{path}{character string specifying the location of the bgx file}
  \item{sep}{separator character (default is tab, \code{"\t"})}
  \item{quote}{character string specifying the quoting characters (disabled by default with \code{quote=""}).  See \code{scan} for further information.}
  \item{header}{logical, TRUE if the bgx file has a header, FALSE otherwise}
  \item{probeStart}{character string, below which the information for
the beads of interest appear. Default value is "\\\[Probes\\\]"}
  \item{controlStart}{character string, below which the information for
the control beads appear. Default value is "\\\[Controls\\\]" }
  \item{...}{further arguments to \code{read.table}}
}

\details{
  The .bgx file is a zip file which contains information about each
  probe on an expression BeadArray.

  To read in the file, you first need to unzip it. To do this, replace
  the .bgx extension with .zip (for example rename
\code{HumanRef-8_V2_0_R0_11223162_A.bgx} as
\code{HumanRef-8_V2_0_R0_11223162_A.zip}) and then unzip this file
(which should leave one file \code{HumanRef-8_V2_0_R0_11223162_A} for
our example).
The unzipped file is tab delimited file and should be readable using
\code{readBGX}. At present this should work for Human and Rat expression
 arrays.  For Mouse arrays, the .bgx has a more complicated structure.
}
\value{
  \code{data.frame} containing information about each bead type (probe
  sequence, ID, control status, etc)
}

\author{Matt Ritchie}

\examples{
#human8bgx = readBGX("HumanRef-8_V2_0_R0_11223162_A", fill=TRUE)
#colnames(human8bgx)
#summary(human8bgx$Status)
#human6bgx = readBGX("HumanWG-6_V2_0_R0_11223189_A", fill=TRUE)
#ratbgx = readBGX("RatRef-12_V1_0_R0_11222119_A", fill=TRUE)
}
\keyword{IO}

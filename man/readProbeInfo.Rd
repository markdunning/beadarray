\name{readProbeInfo}
\alias{readProbeInfo}

\title{Read Probe Information}
\description{
  Given a text file and either \code{BeadLevelList} or
  \code{BeadSummaryList}
  this function reads the specified columns from the file and stores
  them in the genes matrix within the object.
}

\usage{
readProbeInfo(file, beadObject, columns = c("Target", "ProbeID"), header = TRUE, sep = ",")
}

%- maybe also 'usage' for other objects documented here.
\arguments{
\item{file}{File containing information about the genes.}
\item{beadObject}{An object of class \code{BeadLevelList} or
  \code{BeadSummaryList}.}
\item{columns}{Vector specifying the column names of the data that
  should be imported from \code{file}.}
\item{header}{Boolean specifying if \code{file} contains a header row.}
\item{sep}{Character specify the column seperator in \code{file}}
}

\value{Return an object of the same class as \code{beadObject}}

\author{Mike Smith}
\keyword{hplot}

\name{plotRG}
\alias{plotRG}
\title{Plot bead level data: R vs G intensities}
\description{
  Shows R and G intensities fom a BeadLevelList object. The beads can
  either be defined to be all beads with particular ProbeIDs
  or rows BeadIDs
}
\usage{
plotRG(BLData, ProbeIDs=NULL, BeadIDs=NULL, log=TRUE, arrays=1,
                   xlim=c(8,16), ylim=c(8,16), xlab="G intensities",
                   ylab="R intensities",
                   main=arrayNames(BLData)[arrays], smooth=TRUE,
                   cols=NULL, ...)
}

\arguments{
  \item{BLData}{a BeadLevelList object}
  \item{ProbeIDs}{a vector of ProbeIDs to plot}
  \item{BeadIDs}{logical/integer vector specifying which rows of data to
    plot (used if ProbeIDs is NULL)}
  \item{log}{logical, if TRUE, take log base 2 of intensities}
  \item{arrays}{which array/s to plot}
  \item{xlim}{x-axis limits for plot}
  \item{ylim}{y-axis limits for plot}
    \item{xlab}{character string specifying x-axis label}
  \item{ylab}{character string specifying y-axis label}
  \item{main}{main plot title}
  \item{smooth}{logical, whether to smooth the points
    (only used when one array is selected for plotting)}
  \item{cols}{colours to use on the plot}
  \item{\dots}{further graphical parameters to \code{plot}}
}
\details{
  The R and G intensities from selected beads and arrays are plotted. 
}
\value{
Plot is produced on the current graphical device.
}

\examples{
data(BLData)

#plot all beads from array 1
plotRG(BLData, arrays=1)

#plot Beads with ProbeIDs 2 on arrays 1-4
plotRG(BLData, arrays=1:4, ProbeIDs=1000, pch=16)
}

\author{Matt Ritchie}

\keyword{hplot}
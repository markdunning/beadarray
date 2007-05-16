\name{qcBeadLevel}
\alias{qcBeadLevel}
\title{Generate simple diagnostic plots for Illumina bead level data}
\keyword{hplot}
\description{Generate simple diagnostic plots for Illumina bead level data}

\usage{
qcBeadLevel(object, whatToPlot="G", RG=FALSE, log=TRUE, nrow = 100, ncol= 100,
                    colDens=1, colBox=0, html=TRUE,
                    fileName="qcsummary.htm", plotdir=NULL,
                    experimentName=NULL, targets=NULL, ...)
}

\arguments{
  \item{object}{a BeadLevelList object}
  \item{whatToPlot}{character string or vector specifying which
    intensities to plot. Possibilities are \code{"G"}, \code{"Gb"} for
    single channel data and \code{"G"}, \code{"Gb"}, \code{"R"} and
    \code{"Rb"} for two-colour data}
  \item{RG}{if TRUE, plot R vs G intensities per array.  Default value
    is FALSE.  Only useful for two-channel data}
  \item{log}{if TRUE log2 intensities are plotted}
  \item{nrow}{How many rows should the array be divided into (used by imageplot)}
  \item{ncol}{How many columns should the array be divided into (used by
    imageplot)}
  \item{colDens}{colors for density plots (default is 1)}
  \item{colBox}{colors for box plot (default is 0)}
  \item{html}{logical scalar.  If \code{TRUE} an html summary page 
    is generated.  If \code{FALSE}, no summary page is generated.}
  \item{fileName}{name of html summary page.  Default is 'qcsummary.htm'.}
  \item{plotdir}{optional character string specifying the filepath where
    the plots will be saved.  Defaults to current working directory.}
  \item{experimentName}{name to appear on HTML report (default is
    NULL).}
  \item{targets}{data.frame containing smaple information}
  \item{...}{further arguments that can be passed to the plotting functions.}
}

\details{
  This function creates boxplots, smoothed histogram (density) plots and
  imageplots of raw bead level intensity data.
  
  An html page which displays the results, is created when \code{html=TRUE}.  The html page name is specified by the \code{fileName} argument.
}

\author{Matt Ritchie}
\examples{
#data(BLData)
#qcBeadLevel(BLData)
}

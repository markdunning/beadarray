\name{probePairsPlot}
\alias{probePairsPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
probePairsPlot(BLData, array = 1, plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{ ~~Describe \code{BLData} here~~ }
  \item{array}{ ~~Describe \code{array} here~~ }
  \item{plot}{ ~~Describe \code{plot} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (BLData, array = 1, plot = FALSE) 
{
    controls = getControlAnno(BLData)
    t1 <- getArrayData(BLData, what = "G", array = array, log = T)
    t2 <- getArrayData(BLData, what = "ProbeID", array = array, 
        log = T)
    pms = unique(controls$Array_Address[grep("phage_lambda_genome:pm", 
        controls$Reporter_Group_Identifier)])
    mms = unique(controls$Array_Address[grep("phage_lambda_genome:mm2", 
        controls$Reporter_Group_Identifier)])
    pms = pms[pms \%in\% unique(t2)]
    mms = mms[mms \%in\% unique(t2)]
    uniqueSeqs = unique(controls$Probe_Sequence[match(pms, controls$Array_Address)])
    pms = pms[match(uniqueSeqs, controls$Probe_Sequence[match(pms, 
        controls$Array_Address)])]
    uniqueSeqs = unique(controls$Probe_Sequence[match(mms, controls$Array_Address)])
    mms = mms[match(uniqueSeqs, controls$Probe_Sequence[match(mms, 
        controls$Array_Address)])]
    count = 1
    labs = as.character(pms[1])
    pmSeq = controls$Probe_Sequence[match(pms[1], controls$Array_Address)]
    mmSeqs = controls$Probe_Sequence[match(mms, controls$Array_Address)]
    for (j in 1:length(mms)) {
        sequenceDiff = sum(unlist(strsplit(as.character(pmSeq), 
            "")) != unlist(strsplit(as.character(mmSeqs[j]), 
            "")))
        if (sequenceDiff == 2) 
            mm2 = mms[j]
    }
    labs = c(labs, as.character(mm2))
    newinten = t1[t2 == pms[1]]
    ys = newinten
    xs = rep(count, length(newinten))
    cols = rep("red", length(newinten))
    count = count + 1
    newinten = t1[t2 == mm2]
    ys = c(ys, newinten)
    xs = c(xs, rep(count, length(newinten)))
    cols = c(cols, rep("magenta", length(newinten)))
    count = count + 1
    for (i in 2:length(pms)) {
        pmSeq = controls$Probe_Sequence[match(pms[i], controls$Array_Address)]
        mmSeqs = controls$Probe_Sequence[match(mms, controls$Array_Address)]
        for (j in 1:length(mms)) {
            sequenceDiff = sum(unlist(strsplit(as.character(pmSeq), 
                "")) != unlist(strsplit(as.character(mmSeqs[j]), 
                "")))
            if (sequenceDiff == 2) 
                mm2 = mms[j]
        }
        labs = c(labs, as.character(pms[i]))
        newinten = t1[t2 == pms[i]]
        ys = c(ys, newinten)
        xs = c(xs, rep(count, length(newinten)))
        cols = c(cols, rep("red", length(newinten)))
        count = count + 1
        labs = c(labs, as.character(mm2))
        newinten = t1[t2 == mm2]
        ys = c(ys, newinten)
        xs = c(xs, rep(count, length(newinten)))
        cols = c(cols, rep("magenta", length(newinten)))
        count = count + 1
    }
    plot(0, 0, xlim = c(0, count), ylim = c(3, 16), type = "n", 
        axes = F, ylab = "log-intensity", xlab = "", main = "Low Stringency Hyb Control", 
        xaxs = "i")
    axis(2)
    box()
    points(xs + rnorm(length(xs), 0, 0.02), ys, col = cols, pch = 16)
    axis(1, labels = labs, at = 1:(count - 1), las = 2)
    abline(v = seq(0.5, count, by = 2), lty = 2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

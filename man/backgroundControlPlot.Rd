\name{backgroundControlPlot}
\alias{backgroundControlPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
backgroundControlPlot(BLData, array = 1, plot = FALSE)
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
    mrn = controls$Array_Address[controls$Reporter_Group_Name == 
        "negative"]
    mrn = unique(mrn[mrn \%in\% unique(t2)])
    negvals = exprs(createBeadSummaryData(BLData, array = array, 
        probes = mrn))[, 1]
    M = median(negvals)
    MAD = mad(negvals)
    outliers = which(negvals > M + 3 * MAD)
    negvals = negvals[-outliers]
    detect = function(x) 1 - (sum(x > negvals)/(length(negvals)))
    totest = t1[t2 \%in\% mrn]
    output = 100 * length(which(sapply(negvals, detect) < 0.05))/length(negvals)
    output = c(output, mean(negvals[-outliers]), var(negvals[-outliers]))
    output
    if (plot) {
        linmat <- matrix(NA, nrow = length(mrn), ncol = 3)
        t3 <- match(t2, mrn)
        for (i in 1:length(mrn)) {
            linmat[i, ] <- quantile(t1[t2 == mrn[i]], c(0.25, 
                0.5, 0.75), na.rm = TRUE)
        }
        linmat <- linmat[order(linmat[, 2]), ]
        plot(0, 0, type = "n", xlim = c(0, 7.6), ylim = c(min(linmat, 
            na.rm = T), max(linmat, na.rm = T)), main = "negative controls", 
            xlab = "Permuted negative controls", ylab = "IQR of log-intensity", 
            axes = F)
        axis(2)
        box()
        for (i in 1:length(mrn)) {
            lines(c(i/100, i/100), c(linmat[i, 1], linmat[i, 
                3]), col = "grey50")
        }
        lines((1:nrow(linmat))/100, linmat[, 2])
    }
    output
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

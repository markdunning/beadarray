\name{calculateDetection}
\alias{calculateDetection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calculateDetection(BSData)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BSData}{ ~~Describe \code{BSData} here~~ }
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
function (BSData) 
{
    detScores = matrix(nrow = nrow(exprs(BSData)), ncol = ncol(exprs(BSData)))
    negControls = getNegatives(BSData)
    detect = function(x) 1 - (sum(x > negvals)/(length(negvals)))
    for (i in 1:ncol(exprs(BSData))) {
        negvals = negControls[, i]
        M = median(negvals)
        MAD = mad(negvals)
        negvals = negvals[negvals < M + 3 * MAD]
        detScores[, i] = sapply(exprs(BSData)[, i], detect)
    }
    detScores
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

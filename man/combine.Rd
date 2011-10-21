\name{combine}
\docType{methods}
\alias{combine,beadLevelData,beadLevelData-method}
\alias{combine,ExpressionSetIllumina,ExpressionSetIllumina-method}
\title{Combine two objects.}
\description{
	Combine two seperate objects into a single object.
}
\usage{

\S4method{combine}{beadLevelData,beadLevelData}(x, y)

\S4method{combine}{ExpressionSetIllumina,ExpressionSetIllumina}(x,y)

}
\arguments{
\item{x}{An object of class \code{\link[=beadLevelData-class]{beadLevelData}} or \code{\link[=ExpressionSetIllumina-class]{ExpressionSetIllumina}}.}
\item{y}{An object of the same class as \code{x}.}
}

\value{
	Returns an object of the same class as the two inputs.
}

\details{

The \code{combine} function allows two objects of the same class that have been created seperately to be combined into one.  

}

\examples{

if(require(beadarrayExampleData)){

data(exampleBLData)

sectionNames(exampleBLData)

data2 <- combine(exampleBLData, exampleBLData)

sectionNames(data2)

data(exampleSummaryData)

d1 <- exampleSummaryData[,c(1:4)]
d2 <- exampleSummaryData[,c(5:8)]

d3 <- combine(d1, d2)

d3


}

}

\author{Mark Dunning, Mike Smith}


\keyword{methods}

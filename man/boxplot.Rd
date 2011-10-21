\name{boxplot-methods}
\docType{methods}

\alias{boxplot,ExpressionSetIllumina-method}

\title{Boxplots from summary data}

\description{
The standard boxplot function has been extended to work with the \code{ExpressionSetIllumina} class. Moreover, it generates graphics using the \code{ggplot2} package and can incorporate user-defined factors into the plots.

}

\details{
Extra factors can be added to the plots provided they are present in either the \code{phenoData} or \code{featureData} or the object. 
}

\value{
A ggplot object is produced and displayed on screen
}

\examples{
    if(require(beadarrayExampleData)){

	data(exampleSummaryData)
	
	subset <- channel(exampleSummaryData, "G")[,1:8]

	boxplot(subset)

	boxplot(subset, what="nObservations")


	###You can use columns from the featureData in the plots. Here we will use the control-type

	head(fData(subset))

	table(fData(subset)[,"Status"])

	boxplot(subset, probeFactor = "Status")

	###Similarly, we group samples according to colums in phenoData

	pData(subset)


	boxplot(subset, sampleFactor = "SampleFac")

	##Both sample and probe factors can be combined into the same plot

	boxplot(subset, sampleFactor = "SampleFac", probeFactor = "Status")

	##Suppose we have found differentially expressed genes between experimental conditions and want to plot their response. This can be done by first subsetting the ExpressionSetIllumina object and then using the probeFactor and sampleFactor accordingly
      
	if(require(illuminaHumanv3.db)){
	
	   ids <- unlist(mget("ALB", revmap(illuminaHumanv3SYMBOL)))
	   subset2 <- subset[ids,]

	   boxplot(subset2, sampleFactor = "SampleFac")
	   boxplot(subset2, sampleFactor = "SampleFac", probeFactor = "IlluminaID")
	}

    }
}

\author{Mark Dunning}

\keyword{methods}


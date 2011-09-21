plotByProbeFactor <- function(data, probeFactor = "PROBEQUALITY"){

	require("reshape")
	require("ggplot2")

	if(is.character(probeFactor)){

	   if(probeFactor %in% colnames(fData(data))) probeFactor = fData(data)[,probeFactor]
	   else stop("Could not find a factor with name ", probeFactor, " in the featureData for this object\n")


	}

	newDf <- melt(data.frame(exprs(data), ProbeFactor = probeFactor), id.vars ="ProbeFactor")

	#	newDf <- newDf[-which(is.na(newDf$mostGs)),]
	
	p<- ggplot(newDf, aes(x=factor(ProbeFactor),y=value,fill=factor(ProbeFactor))) + geom_boxplot() + facet_wrap(~variable,ncol=min(dim(data)[2], 6))

	p


}

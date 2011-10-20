setMethod("boxplot",
    signature(x = "beadLevelData"),
    function (x, transFun=logGreenChannelTransform,...) 
    {

  	arraynms = sectionNames(x)
  	narrays = length(arraynms)

	inten <- unlist(lapply(1:narrays, function(i) transFun(x, array=i)))
	section <- NULL

	for(i in 1:narrays){

		section <- c(section, rep(sectionNames(x)[i], numBeads(x)[i]))
	}

	df <- data.frame(Value = inten, Section = section)

	p <- ggplot(df, aes(x = factor(Section), y = Value, fill=factor(Section))) + geom_boxplot(outlier.shape=NA)
	p
}
)


setMethod("boxplot",
    signature(x = "ExpressionSetIllumina"),
    function (x, what="exprs",probeFactor = NULL, sampleFactor=NULL, facet = NULL,plot=FALSE,...) 
    {
	
	addedSampleFactor <- addedProbeFactor <- FALSE
	

	data <- melt(assayDataElement(x, what))

	if(!is.null(sampleFactor)){

		if(sampleFactor %in% colnames(pData(x))){
		
		data <- data.frame(data, sampleFactor = pData(x)[match(data[,"X2"], rownames(pData(x))),sampleFactor])
		addedSampleFactor <- TRUE

	

		}
	}	

	if(!is.null(probeFactor)){
		if(probeFactor %in% colnames(fData(x))){
		
		data <- data.frame(data, probeFactor = fData(x)[match(data[,"X1"], rownames(fData(x))),probeFactor])
		addedProbeFactor <- TRUE
		}
	}	


	##Traditional boxplot of all probes on every array vs array.	
	
	if(!(addedSampleFactor) & !(addedProbeFactor)){
	
	p <- ggplot(data, aes(x = factor(X2), y = value, fill=factor(X2))) + geom_boxplot(outlier.shape=NA)

	}

	##Change x axis to vary with the new probe factor	

	else if(addedProbeFactor){
		
		p <- ggplot(data, aes(x = factor(probeFactor), y = value, fill=factor(probeFactor))) + geom_boxplot(outlier.shape=NA)

		if(addedSampleFactor){

			p <- p + facet_wrap(~sampleFactor)
		}

		else p <- p + facet_wrap(~X1)

	}
	

	##change x axis to vary with new sample factor
	else if(addedSamplefactor){

		p <- ggplot(data, aes(x = factor(sampleFactor), y = value, fill=factor(sampleFactor))) + geom_boxplot(outlier.shape=NA)
			
		

	}
	
	p

}
)





calculateDetection = function(BSData, status=fData(BSData)$Status, negativeLabel="negative"){

detScores = matrix(nrow=nrow(exprs(BSData)),ncol=ncol(exprs(BSData)))
 
##Use the fData to find which the negative IDS are 

negInd = which(status == negativeLabel)


detect= function(x) 1 - (sum(x>negvals)/(length(negvals)))


	for(i in 1:ncol(exprs(BSData))){

		negvals = exprs(BSData)[negInd,i]

		M = median(negvals, na.rm=TRUE)
		MAD = mad(negvals, na.rm=TRUE)

		negvals = negvals[negvals < M+3*MAD]
		negvals = negvals[!is.na(negvals)]

		detScores[,i] = sapply(exprs(BSData)[,i], detect)

	}

colnames(detScores) =sampleNames(BSData)
rownames(detScores) = featureNames(BSData)


detScores

}



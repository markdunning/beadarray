
getAnnotation = function(BLData){

BLData@experimentData$annotation

}


setAnnotation = function(BLData, annoName){

data(ExpressionControlData)

if(!annoName %in% names(ExpressionControlData)) stop("Supplied annotation name must be one of:" , paste(names(ExpressionControlData), collapse=" "))

BLData@experimentData$annotation = annoName

BLData

} 



checkPlatform <- function(BLData,verbose=FALSE){

	sigsPath = system.file(package="beadarray", "extdata")	
	load(paste(sigsPath, "/platformSigs.Rda",sep=""))


	ids = getBeadData(BLData, array=1, what="ProbeID")

	rks = sapply(platformSigs,function(x) (sum(ids %in% x$V1)/length(ids))*100)


	if(verbose){
	 cat("Percentage of overlap with IDs on this array and known expression platforms\n")
	 show(rks)
	
	}


	if(all(rks < 90)) warning("Choice of platform may not be accurate. Consider re-running checkPlatform with verbose = TRUE option\n")

	names(sort(rks,decreasing=TRUE)[1])


}




beadStatusVector = function(BLData, array=1, controlProfile = NULL){

	if(is.null(controlProfile)){

		annoName = getAnnotation(BLData)
		
		if(is.null(annoName)) stop("No annotation for this beadLevelData")	

		controlProfile = makeControlProfile(annoName)
	}

			
	tmp = BLData[[array]]

	pIDs = tmp[,1]

	statusVector = rep("regular", length(pIDs))

	controlTypes = unique(controlProfile[,2])

	cIDs = split(controlProfile[,1], controlProfile[,2])

	for(i in 1:length(cIDs)){

		statusVector[which(pIDs %in% cIDs[[i]])] = names(cIDs)[i]
	}

	statusVector		
}


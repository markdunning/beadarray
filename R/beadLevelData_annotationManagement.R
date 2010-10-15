
getAnnotation = function(BLData){

BLData@experimentData$annotation

}


setAnnotation = function(BLData, annoName){

data(ExpressionControlData)

if(!annoName %in% names(ExpressionControlData)) stop("Supplied annotation name must be one of:" , paste(names(ExpressionControlData), collapse=" "))

BLData@experimentData$annotation = annoName

BLData

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

		statusVector[match(cIDs[[i]], pIDs)] = names(cIDs)[i]
	}

	statusVector		
}


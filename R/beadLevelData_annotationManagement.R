
getAnnotation <- function(BLData){

.Deprecated("annotation", package="beadarray")

BLData@experimentData$annotation

}


setAnnotation <- function(BLData, annoName){

.Deprecated("anotation<-", package="beadarray")

data(ExpressionControlData)

if(!annoName %in% names(ExpressionControlData)) stop("Supplied annotation name must be one of:" , paste(names(ExpressionControlData), collapse=" "))

BLData@experimentData$annotation = annoName

BLData

} 


#makeControlProfile <- function(annoName){

#	if(annoName %in% names(ExpressionControlData)){

#		controlProfile = data.frame(ArrayAddress = ExpressionControlData[[as.character(annoName)]][,1], Tag = 	ExpressionControlData[[as.character(annoName)]][,3])
#		controlProfile
#	}

#	else cat("Could not make control profile for annotation:", annoName, "\n")



#}




makeControlProfile <- function(annoName){


    annoLoaded <- require(paste("illumina", annoName, ".db",sep=""), character.only=TRUE)

    if(annoLoaded){
  
    
      mapEnv <-  as.name(paste("illumina", annoName, "REPORTERGROUPNAME",sep=""))

      controlInfo <- unlist(as.list(eval(mapEnv)))

      controlIDs <- names(controlInfo)[controlInfo != ""]

      reporterNames <- controlInfo[controlInfo != ""]

      controlArrayAddress <- unlist(mget(controlIDs, eval(as.name(paste("illumina", annoName, "ARRAYADDRESS",sep="")))))

#      controlProfile <- data.frame(ArrayAddress = controlArrayAddress, Tag = reporterNames)

      repeatedEntries <- which(sapply(reporterNames, function(x) length(grep(",", x, fixed=TRUE))>0))


      if(length(repeatedEntries) > 0){

      newIDs <- NULL
      newTags <- NULL

	for(j in 1:length(repeatedEntries)){


	  tags <- unlist(strsplit(as.character(reporterNames[repeatedEntries[j]]), ","))
	  
	  newTags <- c(newTags, tags)
	  
	  newIDs <- c(newIDs, rep(controlArrayAddress[repeatedEntries[j]],length(tags)))


	}

      controlArrayAddress <- controlArrayAddress[-repeatedEntries]
      reporterNames <- reporterNames[-repeatedEntries]
      
      controlArrayAddress <- c(controlArrayAddress, newIDs)	
      reporterNames <- c(reporterNames, newTags)


      }

      data.frame(ArrayAddress = controlArrayAddress, Tag = reporterNames)
      

  }


}



setMethod("annotation", signature(object = "ExpressionSetIllumina"), function(object) object@annotation)
setMethod("annotation", signature(object = "beadLevelData"), function(object) object@experimentData$annotation)

setReplaceMethod("annotation",
                 signature=signature(
                   object="beadLevelData",
                   value="character"),
                 function(object, value) {
                     object@experimentData$annotation <- value
                     object
                 })

setReplaceMethod("annotation",
                 signature=signature(
                   object="ExpressionSetIllumina",
                   value="character"),
                 function(object, value) {
                     object@annotation <- value
                     object
                 })

 
 


checkPlatform <- function(BLData,verbose=FALSE){

	.Deprecated("suggestAnnotation", package="beadarray")

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


suggestAnnotation <- function(data,verbose=FALSE){

	sigsPath = system.file(package="beadarray", "extdata")	
	load(paste(sigsPath, "/platformSigs.Rda",sep=""))


	ids = getBeadData(data, array=1, what="ProbeID")

	rks = sapply(platformSigs,function(x) (sum(ids %in% x)/length(ids))*100)


	if(verbose){
	 cat("Percentage of overlap with IDs on this array and known expression platforms\n")
	 show(rks)
	
	}


	if(all(rks < 90)) warning("Choice of platform may not be accurate. Consider re-running checkPlatform with verbose = TRUE option\n")

	names(sort(rks,decreasing=TRUE)[1])


}




beadStatusVector <- function(BLData, array=1, controlProfile = NULL){

	if(is.null(controlProfile)){

		annoName <- annotation(BLData)
		
		if(is.null(annoName)) stop("No annotation for this beadLevelData")	

		controlProfile <- makeControlProfile(annoName)
	}

			
	tmp <- BLData[[array]]

	pIDs <- tmp[,1]

	statusVector <- rep("regular", length(pIDs))

	controlTypes <- unique(controlProfile[,2])

	cIDs <- split(controlProfile[,1], controlProfile[,2])

	for(i in 1:length(cIDs)){

		statusVector[which(pIDs %in% cIDs[[i]])] <- names(cIDs)[i]
	}

	statusVector		
}


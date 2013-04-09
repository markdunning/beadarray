###A deprecated function. See suggestAnnotation (below)

checkPlatform <- function(BLData,verbose=FALSE){

    .Deprecated("suggestAnnotation", package="beadarray")

    suggestAnnotation(data = BLData, verbose=verbose)

}


suggestAnnotation <- function(data, verbose = FALSE){

	data(platformSigs)
        
        if(class(data) == "ExpressionSetIllumina") {
            existingAnnotation <- annotation(data)
            if(length(existingAnnotation)) {
                stop("Annotation already set as \"", existingAnnotation, "\"")
            }
            else {
                ids = featureNames(data)
            }
        }
        else if (class(data) == "beadLevelData") {
            ids = getBeadData(data, array=1, what="ProbeID")
        }
        else {
            stop("Argument 'data' is of an unrecognised class")
        }

	rks = sapply(platformSigs,function(x) (sum(ids %in% x)/length(ids))*100)


	if(verbose){
	 cat("Percentage of overlap with IDs on this array and known expression platforms\n")
	 show(rks)
	
	}


	if(all(rks < 90)) warning("Choice of platform may not be accurate. Consider re-running checkPlatform with verbose = TRUE option\n")

	fullname <- tolower(names(sort(rks,decreasing=TRUE)[1]))


	if(length(grep("human", tolower(fullname)) > 0)){

	     vname <-  grep("v", strsplit(as.character(fullname), "")[[1]])

	    shortname <- paste("Humanv", substr(as.character(fullname), vname+1, vname+1),sep="")

	}

	else if(length(grep("mouse", tolower(fullname)) > 0)){

	     vname <-  grep("v", strsplit(as.character(fullname), "")[[1]])

	    shortname <- paste("Mousev", substr(as.character(fullname), vname+1, vname+1),sep="")

	}


	else if(length(grep("rat", tolower(fullname)) > 0)){

	     vname <-  grep("v", strsplit(as.character(fullname), "")[[1]])

	    shortname <- paste("Ratv", substr(as.character(fullname), vname+1, vname+1),sep="")

	}


	shortname


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

##deprecated functions


getAnnotation <- function(BLData){

    .Deprecated("annotation", package="beadarray")

    annotation(BLData)

}


setAnnotation <- function(BLData, annoName){

    .Deprecated("annotation<-", package="beadarray")

    BLData@experimentData$annotation = annoName

    BLData

} 


makeControlProfile <- function(annoName, excludeERCC = TRUE){


    annoLoaded <- require(paste("illumina", annoName, ".db",sep=""), character.only=TRUE)

    if(annoLoaded){
  

      annoPkg <- paste("illumina", annoName, ".db",sep="")

      annoVers <- packageDescription(annoPkg, field="Version")
    
      message(paste("Annotating control probes using package ", annoPkg, " Version:", annoVers, "\n",sep=""))

      mapEnv <-  as.name(paste("illumina", annoName, "REPORTERGROUPNAME",sep=""))

      t <- try(eval(mapEnv),silent=TRUE)

      if(class(t) == "try-error"){
	message(paste("Could not find a REPORTERGROUPNAME mapping in annotation package ", annoPkg,". Perhaps it needs updating?", sep=""))

      }

      else{
      
      controlInfo <- unlist(as.list(eval(mapEnv)))

      #controlIDs <- names(controlInfo)[controlInfo != ""]
      controlIDs <- names(controlInfo)[!is.na(controlInfo)]
      #reporterNames <- controlInfo[controlInfo != ""]
      reporterNames <- controlInfo[!is.na(controlInfo)]

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

      if(excludeERCC){

	if(length(grep("ERCC", reporterNames)) > 0){
	  
	   controlArrayAddress <- controlArrayAddress[-grep("ERCC", reporterNames)]
	   reporterNames <- reporterNames[-grep("ERCC", reporterNames)]
	
	}

      }

      profile <- data.frame(ArrayAddress = controlArrayAddress, Tag = reporterNames)
      



      }


      }

}





 identifyControlBeads <- function(BLData, array=1, controlProfile = NULL){

	if(is.null(controlProfile)){

		annoName <- annotation(BLData)
		
		if(is.null(annoName)) stop("No annotation for this beadLevelData")	

		controlProfile <- makeControlProfile(annoName)
	}


	if(!is.null(controlProfile)){
			
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
    
	else message("Could not identify control beads.\n")
}



beadStatusVector <- function(BLData, array=1, controlProfile = NULL){

	.Deprecated("identifyControlBeads", package="beadarray")

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




addFeatureData <- function(data,toAdd = c("SYMBOL", "PROBEQUALITY", "CODINGZONE", "PROBESEQUENCE"), annotation = NULL){


  ##If we've supplied a character vector, use these to get data from an annotation package 
  if(is(toAdd, "vector")){


   ##if annotation slot is null, assume it is stored with the object

    if(is.null(annotation)){

      ###should use a getAnnotation function when we have one

      annoName <- annotation(data)

    } else {
      annoName <- annotation
    }

    annoLoaded <- require(paste("illumina", annoName, ".db",sep=""), character.only=TRUE)

    if(annoLoaded){
  
      ##should somehow check that the mapping exists!
    
      mapEnv <-  sapply(paste("illumina", annoName, toAdd,sep=""),as.name)

      IDs <- featureNames(data)

      l <- lapply(mapEnv, function(x) mget(IDs, eval(x), ifnotfound = NA))

    

      newAnno <- data.frame(matrix(unlist(l), nrow = length(IDs), byrow=FALSE))
      rownames(newAnno) <- as.character(IDs)
      colnames(newAnno) <- toAdd

    ###merge the myFeatures data frame

    featureData(data) = new("AnnotatedDataFrame", data=data.frame(merge(fData(data), newAnno, by=0,sort=FALSE), row.names=IDs))

     data

    } else {

      stop("Could not load the annotation package ", paste("illumina", annoName, ".db",sep=""))

    }

  }

  else if (is(toAdd, "data.frame")){

       featureData(data) = new("AnnotatedDataFrame", data=data.frame(merge(fData(data), toAdd, by=0,sort=FALSE), row.names=IDs))

	data

  }

  else stop("The toAdd argument must either be a character vector or a data frame\n")
  

}


getPlatformSigs <- function(){

####An internal function to demonstrate how the annotation defintions were generated.
##You will need the libraries lumiHumanIDMapping, lumiMouseIDMapping and lumiRatIDMapping to run this code

#require("lumiHumanIDMapping")

human_conn <- lumiHumanIDMapping_dbconn()

tabs <- dbListTables(human_conn)


platformSigs <- NULL


for(i in 1:length(tabs)){

  x <- tabs[i]

  if("Array_Address_Id" %in% dbListFields(human_conn, as.character(x))){


  platformSigs[[x]] <- as.integer(unlist(dbGetQuery(human_conn, paste("select Array_Address_Id from", as.character(x)))))


  }

  else if ("ProbeId" %in% dbListFields(human_conn, as.character(x))){


  platformSigs[[x]] <- as.integer(unlist(dbGetQuery(human_conn, paste("select ProbeId from", as.character(x)))))


  }

}

#require("lumiMouseIDMapping")

mouse_conn <- lumiMouseIDMapping_dbconn()

tabs <- dbListTables(mouse_conn)



for(i in 1:length(tabs)){

  x <- tabs[i]

  if("Array_Address_Id" %in% dbListFields(mouse_conn, as.character(x))){


  platformSigs[[x]] <- as.integer(unlist(dbGetQuery(mouse_conn, paste("select Array_Address_Id from", as.character(x)))))


  }

  else if ("ProbeId" %in% dbListFields(mouse_conn, as.character(x))){


  platformSigs[[x]] <- as.integer(unlist(dbGetQuery(mouse_conn, paste("select ProbeId from", as.character(x)))))


  }

}

#require("lumiRatIDMapping")

rat_conn <- lumiRatIDMapping_dbconn()

tabs <- dbListTables(rat_conn)



for(i in 1:length(tabs)){

  x <- tabs[i]

  if("Array_Address_Id" %in% dbListFields(rat_conn, as.character(x))){


  platformSigs[[x]] <- as.integer(unlist(dbGetQuery(rat_conn, paste("select Array_Address_Id from", as.character(x)))))


  }

  else if ("ProbeId" %in% dbListFields(rat_conn, as.character(x))){


  platformSigs[[x]] <- as.integer(unlist(dbGetQuery(rat_conn, paste("select ProbeId from", as.character(x)))))


  }

}

platformSigs

}


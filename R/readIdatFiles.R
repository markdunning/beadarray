readIdatFiles <- function(idatFiles = NULL) {
    
    if(is.null(idatFiles)) {
        stop("At least one file must be specified")
    }
    
    ## are there any files with "Red" in the name
    ## if so we should read all the green channel first, then the red.
    twoChannel <- any(grepl("Red", idatFiles))
    if(twoChannel) {
        greenfiles <- idatFiles[grep("Grn", idatFiles)]
        redfiles <- idatFiles[grep("Red", idatFiles)]
        ## make sure they match up
        redfiles <- redfiles[match(strtrim(redfiles, 12), strtrim(greenfiles, 12))]
        idatFiles <- c(greenfiles, redfiles)
    }
    
    ## get the array section names
    sectionNames <- basename(idatFiles)[grep("Grn", basename(idatFiles))]
    sectionNames <- substr(sectionNames, 1, regexpr("_Grn", sectionNames)-1)
    
    BSData <- new("ExpressionSetIllumina")
    
    ## loop through IDAT files, reading one at a time.
    for(file in idatFiles) {
        
        idatData <- readIDAT(file)
        
        ProbeID <- idatData$Quants[,"CodesBinData"]
        
        ## add the expression values, std errors and n observations.  
        ## if this is the first file we stick them straight in, if not we cbind() with existing data
        if( file == idatFiles[1] ) {
            eMat  <- matrix(idatData$Quants[,"MeanBinData"], ncol = 1)
            nObs  <- matrix(idatData$Quants[,"NumGoodBeadsBinData"], ncol = 1)
            varMat <- matrix(idatData$Quants[,"DevBinData"] / sqrt(idatData$Quants[,"NumGoodBeadsBinData"]), ncol = 1)
        }
        else {
            eMat  <- cbind(eMat, idatData$Quants[,"MeanBinData"]) 
            nObs  <- cbind(nObs, idatData$Quants[,"NumGoodBeadsBinData"])  
            varMat  <- cbind(varMat, idatData$Quants[,"DevBinData"] / sqrt(idatData$Quants[,"NumGoodBeadsBinData"]) )    
        }

    }
    assayData(BSData) <- assayDataNew(exprs = eMat, se.exprs = varMat, nObservations = nObs, storage.mode="list")
    
    ## try and work out what platform we've got here
    BSData <- setFeatureData(BSData, ProbeID)
    
    ## set row and column names
    for(index in names(assayData(BSData))) {
        if(length(assayData(BSData)[[ index ]]) > 0) {
            rownames(assayData(BSData)[[ index ]]) <- rownames(pData(featureData(BSData)))
            colnames(assayData(BSData)[[ index ]]) <- sectionNames
        }
    }
    
    ## phenoData
    phenoData(BSData) <- new("AnnotatedDataFrame", data.frame(sectionNames, row.names=sectionNames))
    
    ## set the channel data
    if(twoChannel) {
        channelFac <- c(rep("Green", length(idatFiles)/2, rep("Red", length(idatFiles)/2)))
        channelList <- c("Green", "Red")
    }
    else {
        channelFac <- rep("Green", length(idatFiles))
        channelList <- "Green"
    }
    
    BSData@channelData <- list(channelFac, channelList)
    BSData
}


setFeatureData <- function(BSData, probeIDs) {
    
    annoName <- suggestAnnotation_Vector(probeIDs)
    
    if(!is.null(annoName)) {
 
        annoLoaded <- require(paste("illumina", annoName, ".db",sep=""), character.only=TRUE)
    
        if(annoLoaded){
    
            mapEnv <-  as.name(paste("illumina", annoName, "ARRAYADDRESS",sep=""))
            IlluminaIDs = as.character(unlist(mget(as.character(probeIDs), revmap(eval(mapEnv)),ifnotfound=NA)))
  
            status = rep("Unknown", length(probeIDs)) 
            annoPkg <- paste("illumina", annoName, ".db",sep="")
            annoVers <- packageDescription(annoPkg, fields = "Version")
            
            message(paste("Annotating control probes using package ", annoPkg, " Version:", annoVers, sep=""))
        
            mapEnv <-  as.name(paste("illumina", annoName, "REPORTERGROUPNAME",sep=""))
        
            t <- try(eval(mapEnv),silent=TRUE)
        
            if(class(t) == "try-error"){
                message(paste("Could not find a REPORTERGROUPNAME mapping in annotation package ", annoPkg,". Perhaps it needs updating?", sep=""))    
            }
            else{ 
                status[which(!is.na(IlluminaIDs))] = unlist(mget(IlluminaIDs[which(!is.na(IlluminaIDs))], eval(mapEnv), ifnotfound=NA))	
                status[which(is.na(status))] = "regular"   
            }
            
            ## if we have some ArrayAddressIDs that don't map to an IlluminaID we'll remove the bead-type
            for(index in names(assayData(BSData))) {
                if(length(assayData(BSData)[[ index ]]) > 0) {
                    assayData(BSData)[[ index ]] <- assayData(BSData)[[ index ]][-which(is.na(IlluminaIDs)), ,drop=FALSE ]
                }
            }
            probeIDs <- probeIDs[ -which(is.na(IlluminaIDs)) ]
            status <- status[ -which(is.na(IlluminaIDs)) ]
            IlluminaIDs <- IlluminaIDs[ -which(is.na(IlluminaIDs)) ]
            
            features <- new("AnnotatedDataFrame", data = data.frame(ArrayAddressID = probeIDs, IlluminaID = IlluminaIDs, Status = status, row.names = IlluminaIDs))
        }
    }
    else {
        features <- new("AnnotatedDataFrame", data = data.frame(ArrayAddressID=probeIDs, row.names=probeIDs))
    }
    featureData(BSData) <- features
    annotation(BSData) <- annoName
    return( BSData )
}

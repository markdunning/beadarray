HULK <- function(BLData, array = 1, neighbours = NULL, invasions = 20, useLocs = TRUE, weightName ="wts", transFun = logGreenChannelTransform, outlierFun = illuminaOutlierMethod) {

    message(paste("HULKing section", sectionNames(BLData)[array]))

    if(is.null(neighbours)) {
        message("Calculating Neighbourhood")
        neighbours <- generateNeighbours(BLData, array, useLocs=useLocs)
    }

    data <- transFun(BLData, array)

    return( data - HULKResids(BLData, array = array, transFun = transFun, outlierFun = outlierFun, useLocs, neighbours, invasions, weightName = weightName) )
}

HULKResids <- function(BLData, array, transFun = logGreenChannelTransform, outlierFun = illuminaOutlierMethod, useLocs = TRUE, neighbours = NULL, invasions = 20, weightName = "wts") {

	
    probeIDs = getBeadData(BLData, array = array, what = "ProbeID")
            
    if(weightName %in% colnames(BLData[[array]])){
        weights <- getBeadData(BLData, array = array, what = weightName)
    }
    else {
        weights = rep(1, length(probeIDs))
    }

    data <- transFun(BLData, array)

   ###Remove outliers first	

    outliers = outlierFun(data, probeIDs, wts = weights)
    	
    if(length(outliers)) {
        beadTypeMeans = lapply(split(data[-outliers], probeIDs[-outliers]), mean, na.rm=TRUE)
    }
    else {
        beadTypeMeans = lapply(split(data, probeIDs), mean, na.rm=TRUE)
    }


    ## if an entire bead-type is missing (e.g. all weights 0) then the residuals will all be NA
    residuals = data - unlist(beadTypeMeans)[as.character(probeIDs)]
    ## beads with a residual of 0 are ignored
    residuals[which((is.na(residuals)) | (weights == 0))] = 0
    
    if(is.null(neighbours)) {
        message("Calculating Neighbourhood\n")
        neighbours <- generateNeighbours(BLData, array, useLocs = useLocs)
    }
       
    output <- .C("HULK", as.double(residuals), as.integer(t(neighbours)), as.integer(nrow(neighbours)), as.integer(invasions), results = as.double(residuals))
    output$results
}

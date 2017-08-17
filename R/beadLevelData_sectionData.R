
insertSectionData <- function(BLData, what, data) {


	BLData@sectionData[[what]] <- data.frame(data)
    	return(BLData)


}


getSectionData = function(BLData, what = NULL) {

	if(is.null(what)) BLData@sectionData

	else if(!(what %in% names(BLData@sectionData))) stop("Could not find data ", what, " in sectionData\n, Possible entries", paste(names(BLData@sectionData), collapse=" "))

	else BLData@sectionData[[what]]

}

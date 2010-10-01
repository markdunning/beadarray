##function to add or modify arrayData in a BeadLevelList
insertBeadLevelData <- function(BLData, array = 1, what, data) {
  
    if(array < 1)
        stop("'Array' arguments must be a positive integer");

	secNames = sectionNames(BLData)
	
	#if(array > length(secNames)) stop("Only ", length(secNames), " sections in beadLevelData object. Cannot insert data for array ", array, "\n")

	
	
    arrayName = sectionNames(BLData)[array]

	###Check how many beads for this array and that the data we're trying to insert is the correct length


	##When this function is called by readIllumina, the numBeads will not be set yet
	if(!is.null(numBeads(BLData))){

		if(length(data) != numBeads(BLData)[array]){

		stop("Cannot assign data to this array. Length of data ", length(data), ": Number of beads :", numBeads(BLData)[array], "\n")

		}
	
	}

	BLData@beadData[[arrayName]][[what]] <- new.env()
    	assign(what, data, envir = BLData@beadData[[arrayName]][[what]])
    	return(BLData)


	

}


insertSectionData <- function(BLData, what, data) {


	BLData@sectionData[[what]] <- data.frame(data)
    	return(BLData)


}

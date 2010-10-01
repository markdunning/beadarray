setGeneric("getBeadData", function(BLData, what="Grn", array=1)
   standardGeneric("getBeadData"))

setMethod("getBeadData", "beadLevelData", function(BLData, what="Grn", array=1){
		##Subset to get all data for the array
		tmp = BLData[[array]]

		m = match(what, colnames(tmp))

		if(is.na(m)) stop("Could not find bead data of type ", what, "\n")
	
		tmp[,m]

		}
)





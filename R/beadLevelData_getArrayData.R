setGeneric("getArrayData", function(BLData, what="G", array=1, log=TRUE, method="illumina", n=3, trim=0.05)
   standardGeneric("getArrayData"))

setMethod("getArrayData", "beadLevelData", function(BLData, what="G", array=1,log=TRUE, method="illumina", n=3, trim=0.05){
		##Subset to get all data for the array
		tmp = BLData[[array]]

		m = match(what, colnames(tmp))

		if(is.na(m)) stop("Could not find bead data of type ", what, "\n")
	
		tmp[,m]

		}
)





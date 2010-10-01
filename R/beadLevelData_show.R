## prints a summary of the data contained within a beadLevelData object
## in the same fashion as the original BeadLevelList

setMethod("show", "beadLevelData", function(object) {
   
	cat("Experiment information\n\n")


	show(object@experimentData)	

	ncols = 4
	nrows = 5
	for(i in 1:length(object@sectionData)){

		cat(names(object@sectionData)[i], "\n\n")

		if(is.matrix(object@sectionData[[1]])){
			if(ncol(object@sectionData[[i]]) < ncols){

				show(object@sectionData[[i]][1:nrows,])
				cat(paste("\n...", numBeads(object)[1]-nrows, "more rows of data\n\n"))
			}
			else{
				show(object@sectionData[[i]][1:nrows,1:ncols])
				cat(paste("\n...", numBeads(object)[1]-nrows, "more rows of data\n\n"))
				cat("\n..", ncol(object@sectionData[[i]]) - ncols, "more columns of data\n\n")
			}
		}

		else show(object@sectionData[[i]])

	}

	cat("Array information\n\n")	

	arraynms = sectionNames(object)
	cat(paste("Raw data from section", arraynms[1], "\n\n")) 
	show(object[[1]][1:nrows,])
	cat(paste("\n...", numBeads(object)[1]-nrows, "more rows of data\n\n"))
	cat(paste("... data for", length(arraynms)-1, "more array/s\n\n"))

 })

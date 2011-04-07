## prints a summary of the data contained within a beadLevelData object
## in the same fashion as the original BeadLevelList

setMethod("show", "beadLevelData", function(object) {
   
	cat("\n-----------------\nExperiment information (@experimentData)\n-----------------\n")


	show(object@experimentData)	

	nArrays = length(sectionNames(object))

	ncols = 4
	nrows = 5

    cat("\n-----------------\nPer-section data (@sectionData)\n-----------------\n")        


	for(i in 1:length(object@sectionData)){

		cat(names(object@sectionData)[i], "\n\n")

		if(is.data.frame(object@sectionData[[i]])){
			#if(ncol(object@sectionData[[i]]) < ncols){
            if(ncol(object@sectionData[[i]]) == 1){
				displayrows = nArrays
				show(as.vector(object@sectionData[[i]][1:displayrows,1]))
				cat("\n");
				#if(displayrows > nrows) cat(paste("\n...", nArrays-nrows, "more rows of data\n\n"))
			}
			else{
				displayrows = min(nrows, nArrays)
				show(object@sectionData[[i]][1:displayrows,])

				if(nArrays > displayrows) cat(paste("\n...", nArrays-displayrows, "more rows of data\n\n"))
			}
		}

		else {          
        	show(object@sectionData[[i]])
        }

	}

	cat("\n-----------------\nPer-bead data (@beadData)\n-----------------\n")	
	
	arraynms = sectionNames(object)
	cat(paste("Raw data from section", arraynms[1], "\n\n")) 
	show(object[[1]][1:nrows,])
	cat(paste("\n...", numBeads(object)[1]-nrows, "more rows of data\n\n"))
	cat(paste("... data for", length(arraynms)-1, "more section/s\n\n"))

 })

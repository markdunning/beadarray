## returns the number of beads for all arrays in a beadLevelData object.
## if arrays are specified only numbers for those are returned.

setGeneric("numBeads", function(object, arrays=NULL)
   standardGeneric("numBeads"))

setMethod("numBeads", "beadLevelData", function(object, arrays=NULL) {
    if(is.null(arrays)) 
        object@sectionData$numBeads[,1]
    else 
        object@sectionData$numBeads[arrays,1]
    }
) 

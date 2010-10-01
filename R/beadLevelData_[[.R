## subset a beadLevelData object.
## typing BLData[[1]] returns the beadData relating to the first array

setMethod("[[", 
    signature = c(x = "beadLevelData", i = "ANY", j = "missing"),
    definition = function(x, i, j = "missing") {
        if(length(i)!=1)
            stop("subscript out of bounds (index must have length 1 in '[[')");        
        
        an = sectionNames(x);
        if(i > length(an))
            stop("subscript out of bounds (index greater than the number of arrays)");
        
        out <- NULL
        for(j in names(x@beadData[[an[i]]])) { 
            out <- cbind(out, get(j, x@beadData[[an[i]]][[j]], inherits = FALSE)); 
        }
        colnames(out) <- names(x@beadData[[an[i]]]);
        out
    },
    valueClass = "data.frame")

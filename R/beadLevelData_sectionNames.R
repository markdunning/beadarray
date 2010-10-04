setGeneric("sectionNames", function(object, arrays=NULL)
   standardGeneric("sectionNames"))

setMethod("sectionNames", "beadLevelData", function(object, arrays=NULL) {
   if(is.null(arrays)) 
       as.character(object@sectionData$Targets$sectionName)
   else 
       as.character(object@sectionData$Targets$sectionName[arrays]) 
    }
)
 

setGeneric("metrics", function(object) standardGeneric("metrics"))

setMethod("metrics", signature(object="beadLevelData"), function(object) {
  
    if("Metrics" %in% names(object@sectionData)){

      object@sectionData$Metrics
    }

}

)

setGeneric("p95", function(object,channel) standardGeneric("p95"))

setMethod("p95", signature(object="beadLevelData", channel="character"), function(object,channel="Grn") {
  
    if("Metrics" %in% names(object@sectionData)){

	if(channel == "Grn"){
	  object@sectionData$Metrics$P95Grn
    
        }
	else if(channel == "Red"){
	  object@sectionData$Metrics$P95Red
	} 

	else stop("channel must either Grn or Red")
  
    }

}

)

setGeneric("snr", function(object,channel) standardGeneric("snr"))


setMethod("snr", signature(object="beadLevelData",channel="character"), function(object, channel="Grn") {
  
    if("Metrics" %in% names(object@sectionData)){

	if(channel == "Grn"){
	  object@sectionData$Metrics$P95Grn / object@sectionData$Metrics$P05Grn
    
        }
	else if(channel == "Red"){
	  object@sectionData$Metrics$P95Red / object@sectionData$Metrics$P05Red
	}
#	else stop(channel must either Grn or Red)
  
    }

}

)


setMethod("channelNames", signature(object="ExpressionSetIllumina"), function(object){

	unique(as.character(object@channelData[[1]]))
	}

)

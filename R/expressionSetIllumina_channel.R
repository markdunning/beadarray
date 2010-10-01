setMethod("channel",
    signature(object = "ExpressionSetIllumina"),
    function (object, name, ...) 
    {
        allNames = channelNames(object)

	if(!name %in% allNames) stop("No channel of name ", name, " was not found. Use channelNames function to see what names are valid")

	selArray = which(object@channelData[[1]] == name)	
		
	BSData = new("ExpressionSetIllumina")
	assayData(BSData)=assayDataNew(exprs = exprs(object)[,selArray], se.exprs = se.exprs(object)[,selArray],NoBeads=NoBeads(object)[,selArray],storage.mode="list")
	
	BSData@channelData = object@channelData
	BSData@phenoData = object@phenoData[selArray,]	
	BSData

	}
    
)


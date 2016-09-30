setMethod("channel",
    signature(object = "ExpressionSetIllumina"),
    function (object, name, ...) 
    {
        allNames = channelNames(object)

	if(!name %in% allNames) stop("No channel of name ", name, " was not found. Use channelNames function to see what names are valid")

	selArray = which(object@channelData[[1]] == name)	
		
	BSData = new("ExpressionSetIllumina")
  
	eMat <-exprs(object)[,selArray]
	seMat <- se.exprs(object)[,selArray]
	nobsMat <-  nObservations(object)[,selArray]
	if(!is.null(Detection(object))) detMat <- Detection(object)[,selArray]
	
	newNames <- gsub(paste(name, ":",sep=""), "", colnames(eMat))
	colnames(eMat) <- colnames(seMat) <- colnames(nobsMat) <- newNames
	
	if(!is.null(Detection(object))) colnames(detMat) <- newNames
	
	if(!is.null(Detection(object))){
	  assayData(BSData)=assayDataNew(exprs = eMat, se.exprs = seMat,nObservations=nobsMat,Detection = detMat,storage.mode="list")
	}	else assayData(BSData)=assayDataNew(exprs = eMat, se.exprs = seMat,nObservations=nobsMat,storage.mode="list")
	#assayData(BSData)=assayDataNew(exprs = exprs(object)[,selArray], se.exprs = se.exprs(object)[,selArray],storage.mode="list")

	##Create new channel information

	cData = NULL

	cData[[1]] = object@channelData[[1]][which(object@channelData[[1]] == name)]
	cData[[2]] = object@channelData[[2]][[which(channelNames(object) == name)]]

	BSData@channelData = cData


		


	BSData@phenoData = object@phenoData
	BSData@featureData = object@featureData	
	BSData@QC = object@QC
	annotation(BSData) <- annotation(object)
	BSData
	}
    
)


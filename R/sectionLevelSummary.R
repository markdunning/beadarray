
sectionLevelSummary = function(BLData, channelList, probeIDs=NULL, useMulticore = FALSE, weightNames = "wts", silent=FALSE){

an = sectionNames(BLData)

###If probeIDs not supplied, work out the unique probes in the data

if(is.null(probeIDs)) probeIDs = unique(getBeadData(BLData, what="ProbeID"))

output = vector("list", length(channelList))

###Make template matrices to hold expression values, errors and number of observations


for(cNum in 1:length(channelList)){
	
	template = matrix(nrow=length(probeIDs), ncol=length(an))

	output[[cNum]][["eMat"]] = template

	colnames(output[[cNum]][["eMat"]]) = paste(channelList[[cNum]]@name, an, sep=":")
	rownames(output[[cNum]][["eMat"]]) = probeIDs

	output[[cNum]][["varMat"]]  = template

	colnames(output[[cNum]][["varMat"]]) = paste(channelList[[cNum]]@name, an, sep=":")
	rownames(output[[cNum]][["varMat"]]) = probeIDs


	output[[cNum]][["nObs"]]  = template

	colnames(output[[cNum]][["nObs"]]) = paste(channelList[[cNum]]@name, an, sep=":")
	rownames(output[[cNum]][["nObs"]]) = probeIDs


}

	for(i in 1:length(an)){

		for(ch in 1:length(channelList)){
			chName = channelList[[ch]]@name	
			transFun = channelList[[ch]]@transFun[[1]]
			oFUN = channelList[[ch]]@outlierFun[[1]]
			exprFun = channelList[[ch]]@exprFun[[1]]
			varFun = channelList[[ch]]@varFun[[1]]	
							
			if(!silent) cat("Summarizing ", chName, " channel\n")
					
			if(!silent) cat("Processing Array", i, "\n")
			
			##Make a data-frame of all bead properties on the given array

			tmp = BLData[[i]]	
	
			pidCol = grep("ProbeID", colnames(tmp))	

			###Remove any probes with IDs that are not in 'probeIDs'

			tmp = tmp[which(tmp[,pidCol] %in% probeIDs),]
	

			###Transform to get the values we are interested in summarizing one value per bead
			
			values = transFun(BLData, array=i)		
		
			####Check correct number of values were returned

			if(length(values) != nrow(tmp)) stop("Transformation function did not return correct number of values")

			###Now get vectors of probe IDs and bead weights

			wCol = grep(weightNames, colnames(tmp))
						
			###if no weights were specified, assume all weights are 1			

			if(length(wCol) == 0) wts = rep(1, nrow(tmp))
			
			else wts = tmp[,wCol]	

			pIDs = tmp[,pidCol]
			
		
			###Remove rows that are outliers
			if(!silent) cat("Removing outliers\n")
			oList = oFUN(values, tmp[,pidCol])

			values = values[-oList]
			pIDs = pIDs[-oList]
			wts = wts[-oList]

			###Create list of values, split by ProbeID. Multiply by probe weights

			tmp = split(wts*values, pIDs)
			
			###Find out the mapping between the list and probeIDs

			pMap = match(names(tmp), probeIDs)
	
			if(!silent) cat("Using exprFun\n")
			if(useMulticore) output[[ch]][["eMat"]][pMap,i] = unlist(mclapply(tmp, exprFun))
			else output[[ch]][["eMat"]][pMap,i] = unlist(lapply(tmp, exprFun))
			
			if(!silent) cat("Using varFun\n")		
		
			if(useMulticore) output[[ch]][["varMat"]][pMap,i] = unlist(mclapply(tmp, varFun))			
			
			else output[[ch]][["varMat"]][pMap,i] = unlist(lapply(tmp, varFun))

			if(useMulticore) output[[ch]][["nObs"]][pMap,i] = unlist(mclapply(tmp, length))
			
			else output[[ch]][["nObs"]][pMap,i] = unlist(lapply(tmp, length))

		}

	}


##output

cat("Making  summary object\n")


eMat = output[[1]][["eMat"]]
varMat = output[[1]][["varMat"]]
nObs = output[[1]][["nObs"]]

if(length(output) > 1){

	for(i in 2:length(output)){
		eMat = cbind(eMat, output[[i]][["eMat"]])
		varMat = cbind(varMat, output[[i]][["varMat"]])
		nObs = cbind(nObs, output[[i]][["nObs"]])
	}
}

	channelFac = NULL

	for(i in 1:length(channelList)){

		channelFac = c(channelFac, rep(channelList[[i]]@name, length(an)))
	}


BSData = new("ExpressionSetIllumina")
assayData(BSData)=assayDataNew(exprs = eMat, se.exprs = varMat, NoBeads = nObs,storage.mode="list")

annoName = getAnnotation(BLData)

if(!is.null(annoName)){
	mapPath = system.file(package="beadarray", "extdata")	
	mapTable <- paste(tolower(annoName), "BeadLevelMapping", sep="")
	load(paste(mapPath, "/",mapTable,".rda",sep=""))

	###Must be a way of doing this

	mapTable = eval(as.name(mapTable))	
	
	convertArrayAddressID = function(arrayAddressID, mapTable){
		mapTable[match(arrayAddressID, mapTable[,1]),2]
	}
	
	IlluminaIDs = convertArrayAddressID(probeIDs, mapTable)
	
	isMapped = probeIDs %in% mapTable[,1]	
	


	featureData(BSData) = new("AnnotatedDataFrame", data=data.frame(ProbeID=probeIDs,IlluminaID  = IlluminaIDs, row.names=probeIDs))

}

else featureData(BSData) = new("AnnotatedDataFrame", data=data.frame(ProbeID=probeIDs,row.names=probeIDs))




###We will also set the phenoData slot to be the sectionData from the beadLevelData object


	
BSData@channelData = list(channelFac, channelList)	

BSData

}

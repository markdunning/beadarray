uniqueProbeList = function(BLData){

secNames = sectionNames(BLData)

uIDs = NULL

##Probably only need to look at one sample?

for(i in 1:length(secNames)){

	uIDs = c(uIDs, unique(BLData[[i]][,1]))

}

unique(uIDs)

}

sampleLevelSummary = function(BLData, channelList, probeIDs=NULL, useSampleFac = TRUE, sampleFac= NULL,useMulticore = FALSE, weightNames = "wts"){

arraynms = sectionNames(BLData)



output = vector("list", length(channelList))

###Make template matrices to hold expression values, errors and number of observations


##Get the sample factor from the sectionData slot

sampleFac = getSectionData(BLData, "SampleGroup")[,1]

sList = unique(sampleFac)


##create unique list unless otherwise specified

if(is.null(probeIDs)) {
	cat("Finding list of unique probes in beadLevelData\n")
		
	probeIDs = uniqueProbeList(BLData)
	cat(length(probeIDs), " unique probeIDs found\n")
}


for(cNum in 1:length(channelList)){
	
	template = matrix(nrow=length(probeIDs), ncol=length(sList))

	output[[cNum]][["eMat"]] = template

	colnames(output[[cNum]][["eMat"]]) = paste(channelList[[cNum]]@name, sList, sep=":")
	rownames(output[[cNum]][["eMat"]]) = probeIDs

	output[[cNum]][["varMat"]]  = template

	colnames(output[[cNum]][["varMat"]]) = paste(channelList[[cNum]]@name, sList, sep=":")
	rownames(output[[cNum]][["varMat"]]) = probeIDs


	output[[cNum]][["nObs"]]  = template

	colnames(output[[cNum]][["nObs"]]) = paste(channelList[[cNum]]@name, sList, sep=":")
	rownames(output[[cNum]][["nObs"]]) = probeIDs


}

	currentSamp = sampleFac[1]

	sCount = 1

		

	for(s in 1:length(sList)){

		an = which(sampleFac == sList[s])

		pIDs = wts = NULL

		values = vector("list", length(channelList))

		for(i in an){

			
			tmp = BLData[[i]]	
	
			pidCol = grep("ProbeID", colnames(tmp))	

		
			###Remove any probes with IDs that are not in 'probeIDs'

			tmp = tmp[which(tmp[,pidCol] %in% probeIDs),]
			wCol = grep(weightNames, colnames(tmp))
		
			pIDs = c(pIDs, tmp[,pidCol])
			###If weights were not found, set all weights to 1

			if(length(wCol) == 0) wts = rep(1, length(pIDs))
			
			else wts = c(wts,tmp[,wCol])
			

			for(ch in 1:length(channelList)){

				chName = channelList[[ch]]@name	
				transFun = channelList[[ch]]@transFun[[1]]
				oFUN = channelList[[ch]]@outlierFun[[1]]
				exprFun = channelList[[ch]]@exprFun[[1]]
				varFun = channelList[[ch]]@varFun[[1]]	
							
				cat("Summarizing ", chName, " channel\n")
					
				cat("Processing Array", i, "\n")
			
				###Transform to get the values we are interested in summarizing one value per bead
					
				newVals = transFun(BLData, array=i)					
				
				####Check correct number of values were returned

				if(length(newVals) != nrow(tmp)) stop("Transformation function did not return correct number of values")

				values[[ch]] = c(values[[ch]],newVals)				

			}

		}

		
		for(ch in 1:length(channelList)){
			
			values2 = values[[ch]]
			
			###Remove rows that are outliers

			cat("Removing outliers\n")
			
			###Make sure there are no NA values

			naVals = which(is.na(values2) | is.infinite(values2))

			
			
			values2 = values2[-naVals]
			pIDs2 = pIDs[-naVals]
			wts2 = wts[-naVals]
			
			###Make sure probes and values are ordered according to ProbeID

			pOrder = order(pIDs2)

			pIDs2 = pIDs2[pOrder]
			values2  = values2[pOrder]
			wts2 = wts2[pOrder]

		
			oList = oFUN(values2, pIDs2)

			values2 = values2[-oList]
			pIDs2 = pIDs2[-oList]
			wts2 = wts2[-oList]


			###Create list of values, split by ProbeID. Multiply by probe weights

			tmp = split(wts2*values2, pIDs2)
			
			###Find out the mapping between the list and probeIDs

			pMap = match(names(tmp), probeIDs)
	
			cat("Using exprFun\n")
			output[[ch]][["eMat"]][pMap,s] = unlist(lapply(tmp, exprFun))
			
			cat("Using varFun\n")		
		
			
			output[[ch]][["varMat"]][pMap,s] = unlist(lapply(tmp, varFun))

			
			output[[ch]][["nObs"]][pMap,s] = unlist(lapply(tmp, length))

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

		channelFac = c(channelFac, rep(channelList[[i]]@name, length(sList)))
	}


BSData = new("ExpressionSetIllumina")
assayData(BSData)=assayDataNew(exprs = eMat, se.exprs = varMat, NoBeads = nObs,storage.mode="list")
featureData(BSData) = new("AnnotatedDataFrame", data=data.frame(ProbeID=probeIDs,row.names=probeIDs))

###We will also set the phenoData slot to be the sectionData from the beadLevelData object

BSData@channelData = list(channelFac, channelList)	

BSData

}

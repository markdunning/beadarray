logGreenChannelTransform = function(BLData, array){

	x = getBeadData(BLData, array=array,what="Grn")

    return(log2.na(x))
}

greenChannelTransform = function(BLData, array){

	x = getBeadData(BLData, array=array,what="Grn")
    return(x)

}



logRedChannelTransform = function(BLData,array){

	x = getBeadData(BLData, array=array,what="Red")

    return(log2.na(x))
}


redChannelTransform = function(BLData,array){

	x = getBeadData(BLData, array=array,what="Red")
    return(x)

}


logRatioTransform = function(BLData, array=array){

	x = getBeadData(BLData, array=array,what="Grn")

	y = getBeadData(BLData, array=array,what="Red")	

    return( log2.na(x) - log2.na(y) )	
}




greenChannel <- new("illuminaChannel", logGreenChannelTransform, illuminaOutlierMethod, function(x) mean(x,na.rm=TRUE), function(x) sd(x,na.rm=TRUE),  "G")


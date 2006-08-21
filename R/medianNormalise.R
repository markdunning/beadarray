"medianNormalise" <-
function(BLData){

narrays = ncol(BLData$G)

BLData$G = log2(as.matrix(BLData$G))

med = median(BLData$G,na.rm=TRUE)

for(i in 1:narrays){

BLData$G[,i]  = BLData$G[,i] - median(BLData$G[,i], na.rm=TRUE)

}


BLData$G = BLData$G + med

BLData$normalised = "median"

return(BLData)

}


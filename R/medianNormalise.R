"medianNormalise" <-
function(BLData){

narrays = ncol(BLData$R)

BLData$R = log2(as.matrix(BLData$R))

med = median(BLData$R,na.rm=TRUE)

for(i in 1:narrays){

BLData$R[,i]  = BLData$R[,i] - median(BLData$R[,i], na.rm=TRUE)

}


BLData$R = BLData$R + med

BLData$normalised = "median"

return(BLData)

}


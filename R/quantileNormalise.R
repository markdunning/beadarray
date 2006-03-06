"quantileNormalise" <-
function(BLData, arrays=1:length(BLData$R[1,])){

Rvalues = log2(BLData$R[,arrays])

Rvalues.q = normalizeBetweenArrays(Rvalues, method="quantile")


BLData$R[,arrays] = Rvalues.q

BLData$normalised = "quantile"

BLData

}



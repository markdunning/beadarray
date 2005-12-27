"quantileNormalise" <-
function(BLData){

Rvalues = log2(BLData$R)

Rvalues.q = normalizeBetweenArrays(Rvalues, method="quantile")

BLData.q = BLData

BLData.q$R = Rvalues.q

BLData.q$normalised = "quantile"

BLData.q

}



"quantileNormalise" <-
function(BLData, arrays=1:length(BLData$R[1,])){


  Rvalues = BLData$R[,arrays]
  
Rvalues.q = normalizeBetweenArrays(as.matrix(Rvalues), method="quantile")

BLData$R[,arrays] = Rvalues.q

BLData$normalised = "quantile"

BLData

}



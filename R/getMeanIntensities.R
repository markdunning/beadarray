"getMeanIntensities" <-
function(BSData, probeID){

  if(!class(BSData) == "BeadSummaryList"){
    stop("BeadSummaryList object required!")
  }


len = length(BSData$R[1,])

values = vector(length = len)

BSData$R[BSData$probeID[,1] == probeID,1:len]


}


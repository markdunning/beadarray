"getMeanIntensities" <-
function(BSData, ProbeID){

  if(!class(BSData) == "BeadSummaryList"){
    stop("BeadSummaryList object required!")
  }


len = length(BSData$R[1,])

values = vector(length = len)

BSData$R[BSData$ProbeID == ProbeID,1:len]


}


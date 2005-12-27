"getProbeIndices" <-
function(BLData, probe,array){

all_ids = 1:length(BLData$R[,1])

all_ids[BLData$probeID[,array]==probe]


}


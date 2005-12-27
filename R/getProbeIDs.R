"getProbeIDs" <-
function(BLData){

sort(unique(BLData$probeID[BLData$probeID[,1] > 0,1]))

}


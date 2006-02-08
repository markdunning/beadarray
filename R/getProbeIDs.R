"getProbeIDs" <-
function(BLData){

sort(unique(BLData$ProbeID[BLData$ProbeID[,1] > 0,1]))

}


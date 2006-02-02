"getProbeIndices" <-
function(BLData, probe,array){
which(BLData$probeID[,array] %in% probe)


}


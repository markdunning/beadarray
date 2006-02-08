"getProbeIndices" <-
function(BLData, probe,array){
which(BLData$ProbeID[,array] %in% probe)


}


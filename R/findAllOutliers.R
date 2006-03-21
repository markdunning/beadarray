findAllOutliers <- function(BLData, array, log = FALSE, n = 3, ignoreList=NULL){
  
  probes <- sort(unique(BLData$ProbeID[BLData$ProbeID[,array] > 0,array]))

  intProbeID <- as.integer(sort(BLData$ProbeID[,array]))
  probeIndex <- c(1:length(intProbeID))
  probeIndex <- probeIndex[sort.list(BLData$ProbeID[,array])]
  
  outliers <- findBeadStatus(BLData, array = array,  probes = probes, log = log, n =n, outputValid = FALSE, intProbeID = intProbeID, ignoreList=ignoreList, probeIndex = probeIndex)
  outliers
}


findAllOutliers <- function(BLData, array, log = FALSE, n = 3, ignoreList=NULL){
  
  probes <- unique(BLData$ProbeID[BLData$ProbeID[,array] > 0,array])

  intProbeID <- as.integer(BLData$ProbeID[,array])
  
  outliers <- findBeadStatus(BLData, array = array,  probes = probes, log = log, n =n, outputValid = FALSE, intProbeID = intProbeID, ignoreList=ignoreList)
  outliers
}


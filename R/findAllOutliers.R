findAllOutliers <- function(BLData, array, log = FALSE, n = 3){
  
  probes <- unique(BLData$ProbeID[BLData$ProbeID[,array] > 0,array])

  intProbeID <- as.integer(BLData$ProbeID[,array])
  
  outliers <- findBeadStatus(BLData, array = array,  probes = probes, log = log, n =n, outputValid = FALSE, intProbeID = intProbeID)
  outliers
}

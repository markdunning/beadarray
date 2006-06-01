findAllOutliersSlow <- function(BLData, array, log = FALSE, n = 3, ignoreList=NULL){
  
  probes <- sort(unique(BLData$ProbeID[BLData$ProbeID[,array] > 0,array]))

  intProbeID <- as.integer(sort(BLData$ProbeID[,array]))
  probeIndex <- c(1:length(intProbeID))
  probeIndex <- probeIndex[sort.list(BLData$ProbeID[,array])]
  
  outliers <- findBeadStatus(BLData, array = array,  probes = probes, log = log, n =n, outputValid = FALSE, intProbeID = intProbeID, ignoreList=ignoreList, probeIndex = probeIndex)
  outliers
}

findAllOutliers <- function(BLData, array, log = FALSE, n = 3){

  probes <- sort(unique(BLData$ProbeID[BLData$ProbeID[,array] > 0,array]))

  finten <- BLData$R[,array]
  probeList <- BLData$ProbeID[,array]
  nbeads <- length(BLData$R[,array])

  start = 0

  foo <- .C("findAllOutliers", as.double(finten), binStatus = integer(length = nbeads), as.integer(probeList), as.integer(probes), as.integer(length(probes)), as.integer(nbeads), as.integer(start), PACKAGE = "beadarray")

  which((probeList > 0) & (foo$binStatus == 0))

}

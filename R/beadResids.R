beadResids = function(BLData, which="G", array=1, log=TRUE, n=3) {
  if(which=="G")
    bs = createBeadSummaryData(BLData, log=log, n=n, imagesPerArray=1)@assayData$G[,array,drop=FALSE]
  else # which=="R"
    bs = createBeadSummaryData(BLData, log=log, n=n, imagesPerArray=1)@assayData$R[,array,drop=FALSE]
  ind = match(getArrayData(BLData, which="ProbeID", array=array), rownames(bs))
  resid = getArrayData(BLData, which=which, log=log, array=array)-bs[ind]
}

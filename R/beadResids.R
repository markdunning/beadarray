beadResids = function(BLData, what="G", array=1, log=TRUE, method="illumina", n=3, trim=0.05) {
    bs = createBeadSummaryData(BLData, what=what, log=log, arrays=array, imagesPerArray=1, method=method, n=n, trim=trim)@assayData$exprs[,,drop=FALSE]
    ind = match(getArrayData(BLData, what="ProbeID", array=array), rownames(bs))
    resid = getArrayData(BLData, what=what, log=log, array=array)-bs[ind]
}

#beadResids = function(BLData, which="G", array=1, log=TRUE, n=3) {
#  if(which=="G" & BLData@arrayInfo$channels!="two")
#    bs = createBeadSummaryData(BLData, log=log, n=n, imagesPerArray=1)@assayData$exprs[,array,drop=FALSE]
#  else if(which=="G" & BLData@arrayInfo$channels=="two")
#    bs = createBeadSummaryData(BLData, log=log, n=n, imagesPerArray=1)@assayData$G[,array,drop=FALSE]
#  else # which=="R"
#    bs = createBeadSummaryData(BLData, log=log, n=n, imagesPerArray=1)@assayData$R[,array,drop=FALSE]
#  ind = match(getArrayData(BLData, which="ProbeID", array=array), rownames(bs))
#  resid = getArrayData(BLData, which=which, log=log, array=array)-bs[ind]
#}

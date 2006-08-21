"getProbeIntensities" <-
function(BLData, ProbeIDs, array,log=TRUE){

#Check the object is of class BeadLevelList
  if(class(BLData) != "BeadLevelList"){
    stop("BeadLevelList object required!")
  }
  
if(log){
log2(BLData$G[BLData$ProbeID[,array] %in% ProbeIDs,array])
}
else{
BLData$G[BLData$ProbeID[,array] %in% ProbeIDs,array]
}

}


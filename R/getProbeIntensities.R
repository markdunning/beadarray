"getProbeIntensities" <-
function(BLData, probe, array,log=TRUE){

#Check the object is of class BeadLevelList
  if(class(BLData) != "BeadLevelList"){
    stop("BeadLevelList object required!")
  }
  
if(log){
log2(BLData$R[BLData$probeID[,array]==probe,array])
}
else{
BLData$R[BLData$probeID[,array]==probe,array]
}

}


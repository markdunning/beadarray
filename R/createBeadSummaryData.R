createBeadSummaryData <- function(BLData, log = FALSE, n = 3, arrays=seq(1:length(BLData$G[1,])), imagesPerArray = 2, probes = NULL){

  #Check the object is of class BeadLevelList
  if(class(BLData) != "BeadLevelList"){
    stop("BeadLevelList object required!")
  }

  len = ncol(BLData)

  if(imagesPerArray == 1){
    temp <- BLData[BLData$ProbeID[,1] != 0,1]
  }
  else if(imagesPerArray == 2){
    temp <- rbind(BLData[BLData$ProbeID[,1] != 0,1], BLData[BLData$ProbeID[,2] != 0,2])
  }
  else{
    stop("You can only specify 1 or 2 images per array")
  }
  
  if(is.null(probes)){
    probes = sort(unique(as.vector(temp$ProbeID)))
  }
    probes = probes[probes>0 & !is.na(probes)]
    noprobes = length(probes)

  if(imagesPerArray == 1){
    R = G = Rb = Gb = BeadStDev = NoBeads = nooutliers = matrix(0,nrow = noprobes, ncol=len) }
  else{
     R = G = Rb = Gb = BeadStDev = NoBeads = nooutliers = matrix(0,nrow = noprobes, ncol=(len/2)) }

  i = j = 1
   while(j <= len){
    print(i)
    if(log){
     finten <- log2(temp$G)
     binten <- log2(temp$Gb)
    }
    else {
      finten <- temp$G
      binten <- temp$Gb
    }
     probeIDs <- as.integer(temp$ProbeID)
#    start = (length(which($ProbeID[,i] == 0)))x
     start = 0
     blah <- .C("createBeadSummary",  as.double(finten),  as.double(binten), probeIDs, as.integer(probes), as.integer(noprobes), as.integer(length(temp$G)),
                 fore = double(length = noprobes), back = double(length = noprobes), sd = double(length = noprobes), noBeads = integer(length = noprobes),
                 noOutliers = integer(length = noprobes), nextStart = as.integer(start), PACKAGE = "beadarray")

     G[,i] = blah$fore
     Gb[,i] = blah$back
     BeadStDev[,i] = blah$sd
     NoBeads[,i] = blah$noBeads
     nooutliers[,i] = blah$noOutliers
     j = j+imagesPerArray
     i = i + 1
     rm(probeIDs, blah)
     gc()
     if((imagesPerArray == 1) && (i <= ncol(BLData))){
       temp = BLData[BLData$ProbeID[,i] != 0, i]
     }
     else if((imagesPerArray == 2) && (j < ncol(BLData))){
       temp = rbind(BLData[BLData$ProbeID[,j] != 0,j], BLData[BLData$ProbeID[,j+1] != 0,j+1])
     }
   }
  

  BSData = list()
  BSData$G = G
  BSData$Gb = Gb
  BSData$BeadStDev = BeadStDev
  BSData$NoOutliers = nooutliers
  BSData$Nobeads = NoBeads
  BSData$ProbeID = probes

  class(BSData) = "BeadSummaryList"
  BSData
}

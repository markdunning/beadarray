"findHighestSD" <-
function(BLData, array, limit=1){

#Check the object is of class BeadLevelList
  if(class(BLData) != "BeadLevelList"){
    stop("BeadLevelList object required!")
  }
  
probes = sort(unique(BLData$ProbeID[BLData$ProbeID[,1] > 0,1]))

l = NULL

for(i in 1:length(probes)){

int = getProbeIntensities(BLData, probes[i], array)

if(sd(int, na.rm=TRUE)>limit){

l = c(l, probes[i])

}

}

l

}


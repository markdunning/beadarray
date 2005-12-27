"findLowestCounts" <-
function(BLData, array,  limit=24){

#Check the object is of class BeadLevelList
  if(class(BLData) != "BeadLevelList"){
    stop("BeadLevelList object required!")
  }
  
probes = sort(unique(BLData$probeID[BLData$probeID[,1] > 0,1]))

l = NULL

for(i in 1:length(probes)){

c = getProbeCoords(BLData, probes[i], array)

if (length(c[,1]) < limit){

l = c(l, probes[i] )

}

}

l

}


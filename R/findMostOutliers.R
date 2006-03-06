"findMostOutliers" <-
function(BLData, array, limit=5){

#Check the object is of class BeadLevelList
  if(class(BLData) != "BeadLevelList"){
    stop("BeadLevelList object required!")
  }
  
probes = sort(unique(BLData$ProbeID[BLData$ProbeID[,1] > 0,1]))

l = NULL

for(i in 1:length(probes)){

o = length(findBeadStatus(BLData, probes[i], array))

if(o > limit){

l = c(l, probes[i])

}

}

l

}


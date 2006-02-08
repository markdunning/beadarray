"findHighestChis" <-
function(BLData, array, limit=14){

  #Check the object is of class BeadLevelList
  if(class(BLData) != "BeadLevelList"){
    stop("BeadLevelList object required!")
  }

probes = sort(unique(BLData$ProbeID[BLData$ProbeID[,1] > 0,1]))

l = NULL

for(i in 1:length(probes)){

c = getProbeCoords(BLData, probes[i], array)

r = checkRandomness(BLData, array, c)



if(r > limit){

l = c(l, probes[i])

}



}

l

}


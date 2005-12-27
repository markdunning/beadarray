"findAllOutliers" <-
function(BLData, array, log=FALSE, n=3){

#Check the object is of class BeadLevelList
  if(class(BLData) != "BeadLevelList"){
    stop("BeadLevelList object required!")
  }
  
probes = sort(unique(BLData$probeID[BLData$probeID[,1] > 0,1]))

o=NULL

for(i in 1:length(probes)){

o = c(o, findOutliers(BLData, probes[i], array, log=log, n=n)$outliers)

}

o

}


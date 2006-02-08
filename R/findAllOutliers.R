"findAllOutliers" <-
function(BLData, array, log=FALSE, n=3){

#Check the object is of class BeadLevelList
  if(class(BLData) != "BeadLevelList"){
    stop("BeadLevelList object required!")
  }
  
probes = sort(unique(BLData$ProbeID[BLData$ProbeID[,1] > 0,1]))

o=NULL

for(i in 1:length(probes)){

  print(i)

o = c(o, findOutliers(BLData, probes[i], array, log=log, n=n)$outliers)

}

o

}


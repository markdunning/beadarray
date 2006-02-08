"getProbeCoords" <-
function(BLData, probe, array){

#Finds the co-ordinates of all probes with a particular ID

#Check the object is of class BeadLevelList
  if(class(BLData) != "BeadLevelList"){
    stop("BeadLevelList object required!")
  }
  
#Check that x,y co-ords are available, otherwise stop

  if(is.null(BLData$x)){
    stop("X and Y co-ordinates are not present in this data")
  }
  
arrx <- BLData$x[,array]
arrx<-as.integer(arrx)
arry <- BLData$y[,array]
arry<-as.integer(arry)


#standardise the coordinates to start at 0 


xcoords = (BLData$x[BLData$ProbeID[,array]==probe,array]) - min (arrx) + 1
ycoords = (BLData$y[BLData$ProbeID[,array]==probe,array]) - min (arry) + 1

coords = matrix(nrow=length(xcoords),ncol=2)

coords[,1] = xcoords
coords[,2] = ycoords

coords

}


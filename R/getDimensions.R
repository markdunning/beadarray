"getDimensions" <-
function(BLData, array){

#Check the object is of class BeadLevelList
  if(class(BLData) != "BeadLevelList"){
    stop("BeadLevelList object required!")
  }
  
#Check that x,y co-ords are available

  if(is.null(BLData$x)){
    stop("X and Y co-ordinates are not present in this data")
  }
  
arrx <- BLData$x[,array]
arrx<-as.integer(arrx)
arry <- BLData$y[,array]
arry<-as.integer(arry)

#xmax and ymax are maximum co-ordinates

xmax = max(arrx)-min(arrx)
ymax = max(arry)-min(arry)

dims = list(length=2)

dims$xmax = xmax
dims$ymax = ymax

dims
}


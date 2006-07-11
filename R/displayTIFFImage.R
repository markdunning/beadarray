"displayTIFFImage" <-
function(BLData, array, a=680:720, b=680:720,flip=TRUE,showOutliers=TRUE,locateBeads=FALSE, showUnregistered=FALSE, outliers=NULL){

tif = as.character(BLData$targets$Image1[array])

print(tif)

xt = .Call("readTIFF", tif, PACKAGE = "beadarray")

if(showOutliers & is.null(outliers)){
outliers = findAllOutliers(BLData, array)
}

if(showUnregistered){
u = getUnregisteredBeads(BLData, 1)
}
else u=NULL

if(locateBeads){
showTIFF(xt, BLData, BLData$x[,array], BLData$y[,array], a, b, array=array, outliers=outliers,out=showOutliers,locate="beads",flip=flip, unregs=u)
}
else showTIFF(xt, BLData,BLData$x[,array], BLData$y[,array], a, b, array=array, outliers=outliers,out=showOutliers,flip=flip,unregs=u)


}


"plotBeadIntensities" <-
function(BLData, ProbeIDs, array, n=3,log=FALSE,label=FALSE,main=NULL,ylim=NULL,...){

#Check the object is of class BeadLevelList
  if(class(BLData) != "BeadLevelList"){
    stop("BeadLevelList object required!")
  }


  
inten = getProbeIntensities(BLData, ProbeIDs=ProbeIDs, array, log=log)

me = mean(inten, trim=0.5)
mad = mad(inten)

ids = getProbeIndices(BLData, probe=ProbeIDs, array)

o = findBeadStatus(BLData, ProbeIDs, array)

  
if(!is.null(BLData$x)){

arrx <- BLData$x[,array]
arrx<-as.integer(arrx)
arry <- BLData$y[,array]
arry<-as.integer(arry)

xmax = max(arrx)-min(arrx)
ymax = max(arry)-min(arry)

coords = getProbeCoords(BLData, ProbeIDs, array)

xcentre = round(xmax/2)
ycentre = round(ymax/2)

xcoords = coords[,1]
ycoords = coords[,2]


dist=matrix(nrow=length(xcoords),ncol=1)

dist = sqrt((xcoords-xcentre)^2 + (ycoords - ycentre)^2)

xlab="Distance"
}

else{
dist = 1:length(inten)
xlab="Bead Index"
}

if(is.null(ylim)) ylim=range(0, max(inten)+(n+1)*mad)
if(label) type="n"  
  
plot(dist, inten, xlab=xlab, ylab="Intensity", xlim=range(dist),ylim=ylim,...)


  if(label){
text(dist, inten, seq(1:length(inten)))
}

abline(h=me,col="black")
abline(h=me+n*mad,col="blue")
abline(h=me-n*mad,col="blue")

points(dist[which(ids %in% o)], inten[which(ids %in% o)], col="red", pch=3,...)  
  
}


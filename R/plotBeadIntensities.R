"plotBeadIntensities" <-
function(BLData, probe, array, n=3,log=FALSE,label=FALSE,ylim=NULL,main=NULL,...){

#Check the object is of class BeadLevelList
  if(class(BLData) != "BeadLevelList"){
    stop("BeadLevelList object required!")
  }

inten = getProbeIntensities(BLData, probe, array)

  
if(!is.null(BLData$x)){

arrx <- BLData$x[,array]
arrx<-as.integer(arrx)
arry <- BLData$y[,array]
arry<-as.integer(arry)

xmax = max(arrx)-min(arrx)
ymax = max(arry)-min(arry)

coords = getProbeCoords(BLData, probe, array)

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


o = findOutliers(BLData, probe, array, n=n,log=log)$outliers

dist = dist[!is.na(inten)]
inten = inten[!is.na(inten)]


if(log==TRUE){

if(missing(ylim)){

ylim=range(0,16)

}


if(label){



plot(dist,inten, xlab=xlab, ylab="Intensity", ylim=ylim, xlim=range(dist), type="n", main=main, cex.main=2.1)
text(dist, inten, seq(1:length(inten)))
}


else{
plot(dist,inten, xlab=xlab, ylab="Intensity", ylim=ylim, xlim=range(dist), main=main,cex.main=2.1)
}

abline(h=mean(inten,na.rm=TRUE,trim=0.5),col="black")
abline(h=(mean(inten,na.rm=TRUE, trim=0.5)+n*mad(inten,na.rm=TRUE)),col="blue")
abline(h=(mean(inten,na.rm=TRUE, trim=0.5)-n*mad(inten,na.rm=TRUE)),col="blue")

if(!is.null(BLData$x[o,array])){

xs = BLData$x[o, array] -min(arrx) 
ys = BLData$y[o, array] - min(arry)
ds = sqrt((xs - xcentre)^2 + (ys - ycentre)^2)

points(ds, log2(BLData$R[o, array]), col="red", pch=3)

}
else{
points(dist[which(getProbeIndices(RG, probe, array) %in% o)], log2(BLData$R[o, array]), col="red", pch=3)
}
}



else{

if(missing(ylim)){

maxy = max(mean(2^inten,na.rm=TRUE, trim=0.5) + 4*mad(2^inten, na.rm=TRUE), 2^inten)

ylim=range(0,maxy)

}




if(label){
plot(dist,2^inten, xlab=xlab, ylab="Intensity", ylim=ylim, xlim=range(dist), type="n", main=main,cex.main=2.1)
text(dist,2^inten, seq(1:length(inten)))
}

else{
plot(dist,2^inten, xlab=xlab, ylab="Intensity", ylim=ylim, xlim=range(dist), main=main,cex.main=2.1)
}

abline(h=mean(2^inten,na.rm=TRUE,trim=0.5),col="black")
abline(h=(mean(2^inten,na.rm=TRUE,trim=0.5)+n*mad(2^inten,na.rm=TRUE)),col="blue")
abline(h=(mean(2^inten,na.rm=TRUE,trim=0.5)-n*mad(2^inten,na.rm=TRUE)),col="blue")


if(!is.null(BLData$x[o,array])){

xs = BLData$x[o, array] -min(arrx) 
ys = BLData$y[o, array] - min(arry)
ds = sqrt((xs - xcentre)^2 + (ys - ycentre)^2)

points(ds, log2(BLData$R[o, array]), col="red", pch=3)

}

else{



points(dist[which(getProbeIndices(RG, probe, array) %in% o)], BLData$R[o, array], col="red", pch=3)
}

}







}


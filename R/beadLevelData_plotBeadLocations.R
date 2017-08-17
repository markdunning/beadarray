
plotBeadLocations = function(BLData, ProbeIDs=NULL, BeadIDs=NULL, array=1, SAM=FALSE, xCol = "GrnX", yCol="GrnY", xlab="x-coordinate", ylab="y-coordinate", horizontal = TRUE, main=paste("Bead", ProbeIDs, "locations"),...){


xs.orig = getBeadData(BLData, array=array, what="GrnX")
ys.orig = getBeadData(BLData, array=array, what="GrnY")

##If plotting with the longest edge going along the screen, flip the input coordinates

if(horizontal){

	tmp=ys.orig

	ys = xs.orig

	xs = tmp
	rm(tmp)

}

##If keeping the same orientation as the image, flip the y coordinates so that y=0 is at bottom of screen

else{

	xs = xs.orig
	ys = max(ys.orig) - ys.orig 

}


##xs = max(xs) -xs



##Put the origin at (0,0)

xs = xs - min(xs)
ys = ys - min(ys)

xmax = max(xs)
ymax = max(ys)



plot(1, xlim=range(0:xmax), ylim=range(0:ymax) ,type="n", new=TRUE, axes=FALSE,xlab="",ylab="",xaxs="i",yaxs="i",main=main,...)

if(!(is.null(ProbeIDs))){

	pIDs = getBeadData(BLData, what="ProbeID",array=array)


	xs = xs[which(pIDs %in% ProbeIDs)]
	ys = ys[which(pIDs %in% ProbeIDs)]

}
else{

	xs = xs[BeadIDs]
	ys = ys[BeadIDs]
 
}

##if(SAM) {

##yy = c(ymax/2, 0, 0, ymax/2, ymax, ymax)

##xx = c(0,xmax/4, 0.75*xmax, xmax, 0.75*xmax, xmax/4)

##polygon(xx, yy)

##}

##else polygon(x=c(0,0,xmax,xmax), y=c(0, ymax, ymax,0))


points(xs, ys,...)
box()

}

"plotBeadLocations" <- function(BLData, probeIDs=NULL, beadIDs=NULL, array, label=FALSE,...){

dims = getDimensions(BLData, array)
xmax = dims$xmax
ymax = dims$ymax


if(!is.null(probeIDs)){

xs = ys = NULL

for(i in 1:length(probeIDs)){

xs = c(xs, getProbeCoords(BLData, probeIDs[i], array)[,1])
ys = c(ys, getProbeCoords(BLData, probeIDs[i], array)[,2])


}

coords = matrix(nrow=length(xs), ncol=2)

#standardise the coordinates to start from 0

coords[,1] = xs 
coords[,2] = ys 


plotCoords(xmax, ymax, coords, label=label,...)


}

else{

coords=matrix(nrow=length(beadIDs), ncol=2)

coords[,1] = BLData$x[beadIDs, array] - min(BLData$x[,array])
coords[,2] = BLData$y[beadIDs, array] - min(BLData$y[,array])


plotCoords(xmax, ymax, coords,label=label,...)

}

}
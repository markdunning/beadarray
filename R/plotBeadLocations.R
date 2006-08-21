"plotBeadLocations" <- function(BLData, ProbeIDs=NULL, beadIDs=NULL, array, label=FALSE,...){

dims = getDimensions(BLData, array)
xmax = dims$xmax
ymax = dims$ymax


if(!is.null(ProbeIDs)){

xs = ys = NULL

for(i in 1:length(ProbeIDs)){

xs = c(xs, getProbeCoords(BLData, ProbeIDs[i], array)[,1])
ys = c(ys, getProbeCoords(BLData, ProbeIDs[i], array)[,2])


}

coords = matrix(nrow=length(xs), ncol=2)

#standardise the coordinates to start from 0

coords[,1] = xs 
coords[,2] = ys 


plotCoords(xmax, ymax, coords, label=label,...)


}

else{

coords=matrix(nrow=length(beadIDs), ncol=2)

coords[,1] = BLData$GrnX[beadIDs, array] - min(BLData$GrnX[,array])
coords[,2] = BLData$GrnY[beadIDs, array] - min(BLData$GrnY[,array])


plotCoords(xmax, ymax, coords,label=label,...)

}

}

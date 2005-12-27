"probeCluster" <-
function(BLData, probe, array,main=NULL){

#Function to create a clustering on positions of beads of same type
#BLData is BLData object created by read.maimages
#probe is an identifier for the probe type
#array is the number of the array

coords = getProbeCoords(BLData, probe, array)

plotCluster(coords,main)


}


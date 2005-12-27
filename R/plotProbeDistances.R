"plotProbeDistances" <-
function(BLData, probe, array){

#Plots distances between probes of the same type

coords = getProbeCoords(BLData, probe, array)

plotDistances(coords)


}


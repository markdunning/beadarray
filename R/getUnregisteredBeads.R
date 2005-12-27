"getUnregisteredBeads" <-
function(BLData, array){

l = seq(1:length(BLData$R[,array]))

l[BLData$probeID[,array]==0]

}


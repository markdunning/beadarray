uniqueProbeList = function(BLData){

secNames = sectionNames(BLData)

uIDs = NULL

##Probably only need to look at one sample?

for(i in 1:length(secNames)){

	uIDs = c(uIDs, unique(BLData[[i]][,1]))

}

unique(uIDs)

}

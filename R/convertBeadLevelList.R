## converts an old stlye BeadLevelList into the new beadLevelData

convertBeadLevelList <- function(BeadLevelList) {

    ##check the input is a BeadLevelList
    if(class(BeadLevelList)[1] != "BeadLevelList") 
        stop("This function is designed to convert the old BeadLevelList class into the newer beadLevelData format");
    
    beadLevelData <- new(Class = "beadLevelData");       
    #beadLevelData@history <- appendHistory(beadLevelData, "Created from BeadLevelList object");  
    
    ##convert the beadData 
    sn = BeadLevelList@arrayInfo$arrayNames

    ##manually set the sectionNames for the new object	
    beadLevelData@sectionData$Targets$sectionName = matrix(sn, ncol=1)

    ##loop over each array
    for(i in 1:length(sn)) {

        ##loop over each column of the data
        arrayData <- BeadLevelList@beadData[[ sn[i] ]];
        for(j in 1:ncol(arrayData)) {
            
            ##if the data is full of zeros then just ignore it
            if(length(which(abs(arrayData[,j]) > 0)) == 0) 
                next;
            
            ##change the names of some slots, anything else keep the same
            what <- switch(colnames(arrayData)[j],
                           wts = "Weights",
                           G = "Grn",
                           Gb = "GrnB",
                           R = "Red",
                           Rb = "RedB",
                           colnames(arrayData)[j])                      
            ##add the data
            beadLevelData <- insertBeadData(beadLevelData, array = i, what = what, data = arrayData[,j]);       
        }
    }

    return(beadLevelData);
}

   
        

            

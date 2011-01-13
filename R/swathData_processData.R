processSwathData<-function(twocolour=TRUE, textstring="_perBeadFile.txt", Glocsstring1="-Swath1_Grn.locs", Glocsstring2="-Swath2_Grn.locs", Rlocsstring1="-Swath1_Red.locs", Rlocsstring2="-Swath2_Red.locs", GrnTiffSuffix1 = "-Swath1_Grn.tif", GrnTiffSuffix2 = "-Swath2_Grn.tif", RedTiffSuffix2 = "-Swath2_Red.tif", section_height = 326, section_width = 397, swath_overlap = 81, fullOutput = TRUE, newTextString=".txt", verbose = TRUE){

#writeOutFiles<-function(Swaths,an="arrayname",textstring=".txt",method="Basic"){

    files <- list.files(pattern = textstring)
    arrayNames <- unlist(strsplit(files, textstring))

    for(an in arrayNames) {

        cat(an, "\n");

    ## Read in locs files here to save reading them in in every function
    if(verbose) cat("Reading green locs1... ");
    glocs1 <- BeadDataPackR:::readLocsFile(paste(an, Glocsstring1, sep = ""))
    if(verbose) cat("green locs2... ")
    glocs2 <- BeadDataPackR:::readLocsFile(paste(an, Glocsstring2, sep = ""))
    locslist <- list(glocs1 = glocs1, glocs2 = glocs2);

    if(twocolour){
        if(verbose) cat("red locs1... ")
        rlocs1 <- BeadDataPackR:::readLocsFile(paste(an, Rlocsstring1, sep = ""))
        if(verbose) cat("red locs2... ")
        rlocs2 <- BeadDataPackR:::readLocsFile(paste(an, Rlocsstring2, sep = ""))
        locslist<-list(glocs1=glocs1,glocs2=glocs2,rlocs1=rlocs1,rlocs2=rlocs2)
    }
        
    if(verbose) cat("Done\n")

    # read in text file
    if(twocolour){
        t1 <- matrix(unlist(scan(paste(an, textstring, sep = ""), sep = "\t", what = list(integer(), integer(), numeric(), numeric(), integer(), numeric(), numeric()), skip = 1, quiet = TRUE)), ncol = 7)
    }    
    else {
        t1 <- matrix(unlist(scan(paste(an, textstring, sep = ""), sep = "\t", what = list(integer(), integer(), numeric(), numeric()), skip = 1, quiet = TRUE)), ncol = 4)
    }    

    # Work out which observations are in which swath
    t2 <- assignToImage(t1, an, twocolour = twocolour, locs=locslist, GrnTiffSuffix1 = GrnTiffSuffix1, GrnTiffSuffix2 = GrnTiffSuffix2, verbose = verbose)

    
    Swaths <- genSwaths(t2, an, twocolour = twocolour, GrnTiffSuffix2 = GrnTiffSuffix2, RedTiffSuffix2 = RedTiffSuffix2, section_height = section_height, section_width = section_width, swath_overlap = swath_overlap, locs=locslist,verbose = verbose)

    ## assign each bead to a swath
    writeOutFiles(Swaths, an, newTextString, fullOutput = fullOutput, twocolour = twocolour)

    }

}
    

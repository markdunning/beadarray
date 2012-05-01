processSwathData <- function(inputDir = NULL, outputDir = NULL, twoColour = NULL, textstring="_perBeadFile.txt", segmentHeight = 326, segmentWidth = 397, fullOutput = FALSE, newTextString=".txt", verbose = FALSE){

    ## Hardcode the expected tif and locs file names
    ## if these change we'll need to make them arguments instead
    Glocsstring1 = "-Swath1_Grn.locs";
    Glocsstring2 = "-Swath2_Grn.locs";
    Glocsstring3 = "-Swath3_Grn.locs";  
    Rlocsstring1 = "-Swath1_Red.locs";
    Rlocsstring2 = "-Swath2_Red.locs"; 
    Rlocsstring3 = "-Swath3_Red.locs"; 
    GrnTiffSuffix1 = "-Swath1_Grn.tif";
    GrnTiffSuffix2 = "-Swath2_Grn.tif";
    GrnTiffSuffix3 = "-Swath3_Grn.tif";
    RedTiffSuffix1 = "-Swath1_Red.tif";
    RedTiffSuffix2 = "-Swath2_Red.tif";
    RedTiffSuffix3 = "-Swath3_Red.tif";

    ## if directories weren't specified use the current working directory
    inputDir <- ifelse(is.null(inputDir), getwd(), inputDir);
    outputDir <- ifelse(is.null(outputDir), getwd(), outputDir);
    
    files <- list.files(path = inputDir, pattern = textstring)
    arrayNames <- unlist(strsplit(files, textstring))
    
    ## try and guess the number of colurs by looking for files with "Red" in the name
    if(is.null(twoColour))
        if(numberOfChannels(file.path(inputDir, files[1])) == 2) 
            twoColour = TRUE
        else if(numberOfChannels(file.path(inputDir, files[1])) == 1)
            twoColour = FALSE
        else
            stop("Unexpected file format")

    for(an in arrayNames) {

    cat(an, "\n");
    locslist <- list(Grn = list())
    ## Read in locs files here to save reading them in in every function
    if(verbose) cat("Reading green locs1... ");
    glocs1 <- readLocsFile(file.path(inputDir, paste(an, Glocsstring1, sep = "")))
    locslist$Grn[[1]] <- glocs1;
    if(verbose) cat("green locs2... ")
    glocs2 <- readLocsFile(file.path(inputDir, paste(an, Glocsstring2, sep = "")))
    locslist$Grn[[2]] <- glocs2;
    if(file.exists(file.path(inputDir, paste(an, Glocsstring3, sep = "")))) {
        glocs3 <- readLocsFile(file.path(inputDir, paste(an, Glocsstring3, sep = "")))
        locslist$Grn[[3]] <- glocs3
    }

    if(twoColour){
        locslist$Red <- list()
        if(verbose) cat("red locs1... ")
        locslist$Red[[1]] <- readLocsFile(file.path(inputDir, paste(an, Rlocsstring1, sep = "")))
        if(verbose) cat("red locs2... ")
        locslist$Red[[2]] <- readLocsFile(file.path(inputDir, paste(an, Rlocsstring2, sep = "")))
        if(file.exists(file.path(inputDir, paste(an, Rlocsstring3, sep = "")))) {
            locslist$Red[[3]] <- readLocsFile(file.path(inputDir, paste(an, Rlocsstring3, sep = "")))
        }
    }
        
    if(verbose) cat("Done\n")

    # read in text file
    t1 <- readBeadLevelTextFile(file.path(inputDir, paste(an, textstring, sep = "")))


    # Work out which observations are in which swath
    t2 <- assignToImage(t1, an, inputDir = inputDir, twocolour = twoColour, locslist = locslist, GrnTiffSuffix1 = GrnTiffSuffix1, GrnTiffSuffix2 = GrnTiffSuffix2, verbose = verbose)

    if(length(locslist$Grn) == 3) {
        Swaths <- genThreeSwaths(t2, sectionName = an, inputDir = inputDir, twocolour = twoColour, locsList = locslist, GrnTiffSuffix2 = GrnTiffSuffix2, RedTiffSuffix2 = RedTiffSuffix2, segmentHeight = segmentHeight, verbose = verbose)
    } 
    else {
        Swaths <- genTwoSwaths(t2, sectionName = an, inputDir = inputDir, twocolour = twoColour, locsList = locslist, GrnTiffSuffix2 = GrnTiffSuffix2, RedTiffSuffix2 = RedTiffSuffix2, segmentHeight = segmentHeight, segmentWidth = segmentWidth, verbose = verbose)
    }
    

    ## assign each bead to a swath
    writeOutFiles(Swaths, an, newTextString, fullOutput = fullOutput, twocolour = twoColour, outputDir = outputDir)

    }

}
    

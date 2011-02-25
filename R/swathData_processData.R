processSwathData <- function(inputDir = NULL, outputDir = NULL, twoColour = NULL, textstring="_perBeadFile.txt", segmentHeight = 326, segmentWidth = 397, swathOverlap = 81, fullOutput = FALSE, newTextString=".txt", verbose = FALSE){

    ## Hardcode the expected tif and locs file names
    ## if these change we'll need to make them arguments instead
    Glocsstring1 = "-Swath1_Grn.locs";
    Glocsstring2 = "-Swath2_Grn.locs"; 
    Rlocsstring1 = "-Swath1_Red.locs";
    Rlocsstring2 = "-Swath2_Red.locs"; 
    GrnTiffSuffix1 = "-Swath1_Grn.tif";
    GrnTiffSuffix2 = "-Swath2_Grn.tif";
    RedTiffSuffix2 = "-Swath2_Red.tif";

    ## if directories weren't specified use the current working directory
    inputDir <- ifelse(is.null(inputDir), getwd(), inputDir);
    outputDir <- ifelse(is.null(outputDir), getwd(), outputDir);
    
    ## try and guess the number of colurs by looking for files with "Red" in the name
    if(is.null(twoColour))
        twoColour <- ifelse(any(grepl("Red", list.files(path = inputDir))), TRUE, FALSE);

    files <- list.files(path = inputDir, pattern = textstring)
    arrayNames <- unlist(strsplit(files, textstring))

    for(an in arrayNames) {

    cat(an, "\n");

    ## Read in locs files here to save reading them in in every function
    if(verbose) cat("Reading green locs1... ");
    glocs1 <- BeadDataPackR:::readLocsFile(file.path(inputDir, paste(an, Glocsstring1, sep = "")))
    if(verbose) cat("green locs2... ")
    glocs2 <- BeadDataPackR:::readLocsFile(file.path(inputDir, paste(an, Glocsstring2, sep = "")))
    locslist <- list(glocs1 = glocs1, glocs2 = glocs2);

    if(twoColour){
        if(verbose) cat("red locs1... ")
        rlocs1 <- BeadDataPackR:::readLocsFile(file.path(inputDir, paste(an, Rlocsstring1, sep = "")))
        if(verbose) cat("red locs2... ")
        rlocs2 <- BeadDataPackR:::readLocsFile(file.path(inputDir, paste(an, Rlocsstring2, sep = "")))
        locslist<-list(glocs1=glocs1, glocs2=glocs2, rlocs1=rlocs1, rlocs2=rlocs2)
    }
        
    if(verbose) cat("Done\n")

    # read in text file
    if(twoColour){
        t1 <- matrix(unlist(scan(file.path(inputDir, paste(an, textstring, sep = "")), sep = "\t", what = list(integer(), integer(), numeric(), numeric(), integer(), numeric(), numeric()), skip = 1, quiet = TRUE)), ncol = 7)
    }    
    else {
        t1 <- matrix(unlist(scan(file.path(inputDir, paste(an, textstring, sep = "")), sep = "\t", what = list(integer(), integer(), numeric(), numeric()), skip = 1, quiet = TRUE)), ncol = 4)
    }    

    # Work out which observations are in which swath
    t2 <- assignToImage(t1, an, inputDir = inputDir, twocolour = twoColour, locs=locslist, GrnTiffSuffix1 = GrnTiffSuffix1, GrnTiffSuffix2 = GrnTiffSuffix2, verbose = verbose)

    
    Swaths <- genSwaths(t2, an, twocolour = twoColour, GrnTiffSuffix2 = GrnTiffSuffix2, RedTiffSuffix2 = RedTiffSuffix2, section_height = segmentHeight, section_width = segmentWidth, swath_overlap = swathOverlap, locs = locslist, verbose = verbose)

    ## assign each bead to a swath
    writeOutFiles(Swaths, an, newTextString, fullOutput = fullOutput, twocolour = twoColour, outputDir = outputDir)

    }

}
    

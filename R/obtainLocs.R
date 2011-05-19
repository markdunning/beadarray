## takes the name of a .locs file and returns a two column matrix of coordinates

obtainLocs <- function(fileName, filePath) {
    
    if(grepl(".bab", fileName)) {
        allLocs <- BeadDataPackR:::extractLocsFile(inputFile = fileName, path = filePath);
        ## do we want the locs from the red or green image
        if(grepl("Red", fileName))
            locs <- allLocs[,3:4]
        else 
            locs <- allLocs[,1:2]
    }
    else {
        locs <- BeadDataPackR:::readLocsFile(file.path(filePath, fileName));
    }
    return(locs)
}
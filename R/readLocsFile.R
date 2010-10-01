## takes the name of a .locs file and returns a two column matrix of coordinates

readLocsFile <- function(fileName) {
    
    con <- file(fileName, "rb")
    
    ##read the first two bytes
    readBin(con, integer(), n = 2, size = 4)
    ##3rd bytes tells you how many probes there are
    nprobes <- readBin(con, integer(), n = 1, size = 4)
    
    coords <- matrix(ncol = 2, nrow = nprobes)
    colnames(coords) <- c("X", "Y")
    ##read in the whole file
    tmp <- readBin(con, double(), n = 2*nprobes, size = 4)
    ##store the x and y coords in the two columns
    coords[,1] <- tmp[seq(1, 2*nprobes, 2)]
    coords[,2] <- tmp[seq(2, 2*nprobes, 2)]

    close(con)
    return(coords)
}
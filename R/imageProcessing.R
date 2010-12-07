## Wrapper functions to call the illumina image processing routines.

illuminaForeground <- function(pixelMatrix, beadCoords) {
    ## if a single coordinate pair is passed, it needs to be coerced into a matrix
    beadCoords <- matrix(beadCoords, ncol = 2)
    fg <- .Call("illuminaForeground", pixelMatrix, beadCoords, PACKAGE = "beadarray");
    return(fg);
}

illuminaBackground <- function(pixelMatrix, beadCoords) {
    ## if a single coordinate pair is passed, it needs to be coerced into a matrix
    beadCoords <- matrix(beadCoords, ncol = 2)
    bg <- .Call("illuminaBackground", pixelMatrix, beadCoords, PACKAGE = "beadarray");
    return(bg);
}

illuminaSharpen <- function(pixelMatrix) {
    sh <- .Call("illuminaSharpen", pixelMatrix, PACKAGE = "beadarray");
    return(sh);
}

medianBackground <- function(pixelMatrix, beadCoords) {
    ## if a single coordinate pair is passed, it needs to be coerced into a matrix
    beadCoords <- matrix(beadCoords, ncol = 2)
    bg <- .Call("medianBackground", pixelMatrix, beadCoords, PACKAGE = "beadarray");
    return(bg);
}

## Wrapper functions to call the illumina image processing routines.

illuminaForeground <- function(pixelMatrix, beadCoords) {
    fg <- .Call("foreground", pixelMatrix, beadCoords, PACKAGE = "beadarray");
    return(fg);
}

illuminaBackground <- function(pixelMatrix, beadCoords) {
    bg <- .Call("background", pixelMatrix, beadCoords, PACKAGE = "beadarray");
    return(bg);
}

illuminaSharpen <- function(pixelMatrix) {
    sh <- .Call("sharpen", pixelMatrix, PACKAGE = "beadarray");
    return(sh);
}

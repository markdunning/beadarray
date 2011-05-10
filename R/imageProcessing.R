## Wrapper functions to call the illumina image processing routines.

illuminaForeground <- function(pixelMatrix, beadCoords) {
    ## if a single coordinate pair is passed, it needs to be coerced into a matrix
    if(length(as.matrix(beadCoords)) == 2)
        beadCoords <- matrix(beadCoords, ncol = 2)
    else
        beadCoords <- as.matrix(beadCoords)

    ## we need to see if this is an integer or numeric matrix
    integerBool <- ifelse(class(pixelMatrix[1,1]) == "integer", 1L, 0L)

    fg <- .Call("illuminaForeground", pixelMatrix, beadCoords, integerBool, PACKAGE = "beadarray");
    return(fg);
}

illuminaBackground <- function(pixelMatrix, beadCoords) {
    ## if a single coordinate pair is passed, it needs to be coerced into a matrix
    if(length(as.matrix(beadCoords)) == 2)
        beadCoords <- matrix(beadCoords, ncol = 2)
    else
        beadCoords <- as.matrix(beadCoords)

    ## we need to see if this is an integer or numeric matrix
    integerBool <- ifelse(class(pixelMatrix[1,1]) == "integer", 1L, 0L)

    bg <- .Call("illuminaBackground", pixelMatrix, beadCoords, integerBool, PACKAGE = "beadarray");
    return(bg);
}

illuminaSharpen <- function(pixelMatrix) {
    sh <- .Call("illuminaSharpen", pixelMatrix, PACKAGE = "beadarray");
    return(sh);
}

medianBackground <- function(pixelMatrix, beadCoords) {
    ## if a single coordinate pair is passed, it needs to be coerced into a matrix
    if(length(as.matrix(beadCoords)) == 2)
        beadCoords <- matrix(beadCoords, ncol = 2)
    else
        beadCoords <- as.matrix(beadCoords)

    bg <- .Call("medianBackground", pixelMatrix, beadCoords, PACKAGE = "beadarray");
    return(bg);
}

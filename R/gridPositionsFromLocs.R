## This function matches the values in the locs file to those in the txt, removing those
## which weren't decoded.  The index of a bead within the locs file is then used 
## to find it's position within the array's grid structure.

## It takes the pre-read coords from the txt file, the names of the locs file,
## the dimensions of a segment, hopeully taken from the SDF and the size of the gap to insert
## between segments

neighboursFromLocs <- function(txtFileCoords, locsName, locsPath, nrowPerSegment, nSegs) {
    
    ## read the locs file 
    ## we can move the reading of the locs file outside this function if required
    locsFileCoords <- obtainLocs(locsName, locsPath);
    
#     grid <- .Call("locsIndicesToGrid", as.integer(1:(nrowPerSegment * ncolPerSegment)), as.integer(c(nrowPerSegment, ncolPerSegment, 1)), PACKAGE = "beadarray");
#     
#     tmp <- cbind(1:nrow(grid), grid)
#     tmp2 <- cbind(tmp[,1], t(apply(tmp, 1, FUN = function(x, nrow, ncol) {
#                     tmp = 2 * 0^(x[3] %% 2)
#                     tmp2 <- c(x[1] - nrow, x[1] - nrow - 1 + tmp, x[1] - 1, x[1] + 1, x[1] + nrow, x[1] + nrow - 1 + tmp);
#                     
#                     idx <- NULL
#                     if(any(tmp2 < 1)) 
#                             idx <- c(idx, which(tmp2 < 1))
#                     if(any(tmp2 > (nrow * ncol)) )
#                             idx <- c(idx, which(tmp2 > (nrow * ncol)))
#                     if( (x[2] == 1) | (x[2] == 2) ) {
#                         if(any( (tmp2 %% nrow) == 0 ) )
#                             idx <- c(idx, which((tmp2 %% nrow) == 0 ))
#                     }
#                     else if( (x[2] == 2*nrow-1) | (x[2] == 2*nrow ) ) {
#                         if(any( (tmp2 %% nrow) == 1 ) )
#                             idx <- c(idx, which((tmp2 %% nrow) == 1 ))
#                     }
#                     
#                     if(!is.null(idx))
#                             tmp2[idx] <- NA
#                     return(tmp2)                           
#     }, nrow = nrowPerSegment, ncol = ncolPerSegment)) )

    ncolPerSegment <- nrow(locsFileCoords) / (nrowPerSegment * nSegs);

    segmentNeighbours <- .Call("neighboursFromLocs", as.integer(c(nrowPerSegment, ncolPerSegment)), PACKAGE = "beadarray")

    #nSegs <- nrow(locsFileCoords) / (nrowPerSegment * ncolPerSegment)

    locsNeighbours <- matrix(nrow = nSegs * nrowPerSegment * ncolPerSegment, ncol = 7)
    for(i in 1:nSegs) {
        locsNeighbours[(1 + (i-1) * (nrowPerSegment * ncolPerSegment) ):(i * (nrowPerSegment * ncolPerSegment) ),] <- segmentNeighbours + ( (i-1) * (nrowPerSegment * ncolPerSegment) )
    }

    ## round the values so they match those in the txt file
    locsFileCoords <- matrix(.Call("roundLocsFileValues", locsFileCoords, PACKAGE = "BeadDataPackR"), ncol = 2);
    ## add an indexing column to both sets of coords so we can reorder later
    locsFileCoords <- cbind(1:nrow(locsFileCoords), locsFileCoords);
    txtFileCoords <- cbind(1:nrow(txtFileCoords), txtFileCoords);
    
    ## create keys to match the entries
    ## multiplication might produce the same value twice, some maybe it isn't ideal!!!!
    locsKeyMult <- locsFileCoords[,2] * locsFileCoords[,3];
    txtKeyMult <- txtFileCoords[,2] * txtFileCoords[,3];
    if(any(duplicated(locsKeyMult))) {
        locsKeyDiv <- locsFileCoords[,2] / locsFileCoords[,3]
        txtKeyDiv <- txtFileCoords[,2] / txtFileCoords[,3]
        idx <- which( (locsKeyDiv %in% txtKeyDiv) & (locsKeyMult %in% txtKeyMult) )
    }
    else {
        idx <- which( locsKeyMult %in% txtKeyMult );   
    }

    ## if there's still some duplicated use string concatonation. It's really slow!!
    if(length(idx) != nrow(txtFileCoords)) {
        locsKey <- paste(locsFileCoords[,2], locsFileCoords[,3]);
        txtKey <- paste(txtFileCoords[,2], txtFileCoords[,3]);
        idx <- which(locsKey %in% txtKey);
    }

    ## find which are in both file and remove those that aren't from the locsCoords
    nonDecodedCoords <- locsFileCoords[-idx,];
    locsFileCoords <- locsFileCoords[idx,];
    txtNeighbours <- locsNeighbours[idx,];    
    
    ## make both sets of coordinates in the same order
    ## this is the slowest step by a long way, can it be avoided??
    txtFileCoords <- txtFileCoords[order(txtFileCoords[,2], txtFileCoords[,3]),];
    locsFileCoords <- locsFileCoords[order(locsFileCoords[,2], locsFileCoords[,3]),];
            
    hash <- rbind(cbind(locsFileCoords[,1], txtFileCoords[,1]), cbind(nonDecodedCoords[,1], rep(0, nrow(nonDecodedCoords)) ));
    hash <- hash[order(hash[,1]),]

    neighbours <- .Call("hashLocsToTxtIndices", txtNeighbours, hash[,2], PACKAGE = "beadarray")                                
                                    
    neighbours <- neighbours[order(neighbours[,1]),2:7];
    
    return(neighbours);                                    
}
    
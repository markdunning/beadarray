## functions to convert between position in the linear locs file
## and pair of coordinates in the grid

locsIndicesToGrid <- function(x, nrow, ncol) {
    idx <- x - 1;
    col <- abs( ( ( idx %% (ncol * nrow) ) %/% nrow)) + 1;
    row <- 2 * (idx %% nrow) + 1;    
    if(col %% 2 == 0) 
            row <- row + 1;    
    return(c(col, row));
}


gridToLocsIndices <- function(x, nrow, ncol) {
    col <- x[1]; row <- x[2];
    if(col %% 2 == 1)
        row <- row + 1;
    row <- row %/% 2;
    idx <- ((col - 1) * nrow) + row;
    return(idx)
}

gridShift <- function(seg, xShift, yShift, nRows, nCols) {

    if((xShift %% 2) == (yShift %% 2)) {
        originalGrid <- t(sapply(seg[,5], locsIndicesToGrid, nrow = nRows, ncol = nCols));
        newGrid <- cbind(originalGrid[,1] + xShift, originalGrid[,2] + yShift)
        newIndex <- apply(newGrid, 1, gridToLocsIndices, nRows, nCols)
                
        ## fix those outside the grid to zero
        newIndex[which( (newGrid[,1] < 1) | (newGrid[,2] < 1) )] <- 0
        newIndex[which( (newGrid[,1] > nCols) | (newGrid[,2] > (nRows * 2)) )] <- 0
        
        tmpSeg <- seg
        tmpSeg[which(newIndex == 0),1] <- 0
        tmpSeg[which(newIndex != 0),1] <- seg[newIndex,1]
        
        if(length(which(tmpSeg[,1] == 0))) {
            tmpSeg <- tmpSeg[-which(tmpSeg[,1] == 0),]
        }
        if(length(which(tmpSeg[,2] == 0))) {
            tmpSeg <- tmpSeg[-which(tmpSeg[,2] == 0),]
        }
        return(tmpSeg)  
    }
    else { return(NULL) }
}

 ## move the grid of bead intensities and calculate the mean 
 ## within bead-type variance
testGridShift <- function(seg, xShift, yShift, nRows, nCols) {  

    ## only perform a shift if we are actually doing one
    if((xShift != 0) && (yShift != 0)) {
        tmpSeg <- gridShift(seg, xShift, yShift, nRows, nCols)
    }
    else {
        tmpSeg <- seg
    }
    if(!is.null(tmpSeg))    {
        tmpSeg[which(tmpSeg[,2] <= 0),2] <- 0.001
        s <- split(log2(tmpSeg[,2]), tmpSeg[,1])
        v <- unlist(lapply(s, var, na.rm = TRUE))
        return(mean(v, na.rm = T));
    }
    else { return(NA) }
}



checkRegistration <- function(BLData, array = 1) {
        
    sdfFileName <- file.path(BLData@sectionData$Targets$directory[1], list.files(as.character(BLData@sectionData$Targets$directory[1]), pattern = ".sdf")[1]);
    if(file.exists(sdfFileName)) {
        sdf <- simpleXMLparse(readLines(sdfFileName, warn = FALSE))
    }
    else {
        stop("sdf file cannot be located.\nAborting registration check\n");
    }
    
    nSegs <- as.integer(sdf$RegistrationParameters$SizeBlockY[[1]]);
    nRows <- as.integer(sdf$RegistrationParameters$SizeGridX[[1]]);
    nCols <- as.integer(sdf$RegistrationParameters$SizeGridY[[1]]);
    beadsPerSeg <- nRows * nCols;
    
    res <- list();
    
    for(i in array) {
        
        ## make sure they've specified an array that exists
        if( (i < 0) || (i > nrow(BLData@sectionData$Targets)) ) {
            message(paste("Cannot find array", i, "\n"))
            next;
        }
        
        sectionName <- as.character(BLData@sectionData$Targets[i,"sectionName"]);
        cat(sectionName, "\n");
        res[[sectionName]] <- list("Grn" = NULL);
        
        locs <- readLocsFile(file.path(BLData@sectionData$Targets[i,1], paste(sectionName, "_Grn.locs", sep = "")));
        
        ## remove reliance on BDPR later
        comb <- BeadDataPackR:::combineFiles(BLData[[i]][,c("ProbeID", "Grn", "GrnX", "GrnY")], locs)
        comb <- comb[order(comb[,5]),]
        
        for(j in 1:nSegs) {
            seg <- comb[(((j-1) * beadsPerSeg) + 1):(j * beadsPerSeg),];
            meanVar <- testGridShift(seg, 0, 0, nRows = nRows, nCols = nCols);
            res[[sectionName]][["Grn"]] <- c(res[[sectionName]][["Grn"]], meanVar);
        }
        
        ## detect red channel and check registration if appropriate
        if("Red" %in% colnames(BLData[[1]])) {
            res[[sectionName]][["Red"]] <- NULL
            
            locs <- readLocsFile(file.path(BLData@sectionData$Targets[i,1], paste(sectionName, "_Red.locs", sep = "")));
            ## remove reliance on BDPR later
            comb <- BeadDataPackR:::combineFiles(BLData[[i]][,c("ProbeID", "Red", "RedX", "RedY")], locs)
            comb <- comb[order(comb[,5]),]
            
            for(j in 1:nSegs) {
                seg <- comb[(((j-1) * beadsPerSeg) + 1):(j * beadsPerSeg),];
                meanVar <- testGridShift(seg, 0, 0, nRows = nRows, nCols = nCols);
                res[[sectionName]][["Red"]] <- c(res[[sectionName]][["Red"]], meanVar);
            }
        }
        
    }
    return(res);
}

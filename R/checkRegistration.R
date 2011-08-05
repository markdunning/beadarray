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
        #tmpSeg[which(tmpSeg[,2] <= 0),2] <- 0.001
        s <- split(log2.na(tmpSeg[,2]), tmpSeg[,1])
        v <- unlist(lapply(s, var, na.rm = TRUE))
        
        s2 <- split(log2.na(tmpSeg[sample(1:nrow(tmpSeg), size = nrow(tmpSeg)),2]), tmpSeg[,1])
        v2 <- unlist(lapply(s2, var, na.rm = TRUE))
        #return(mean(v, na.rm = T));
        return(v2 - v)
    }
    else { return(NA) }
}

## calculate the within bead-type variances and 
## compare with randomised IDs
scoreRegistration <- function(seg) {  

    if(!is.null(seg))    {
        s <- split(log2.na(seg[,2]), seg[,1])
        v <- unlist(lapply(s, var, na.rm = TRUE))
        
        s2 <- split(log2.na(seg[sample(1:nrow(seg), size = nrow(seg)),2]), seg[,1])
        v2 <- unlist(lapply(s2, var, na.rm = TRUE))
        return(v2 - v)
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
    coords <- list();
    corners <- list();
    p95 <- NULL;
    tiffs <- NULL;
    metrics <- NULL;
    
    for(i in array) {
        
        ## make sure they've specified an array that exists
        if( (i < 0) || (i > nrow(BLData@sectionData$Targets)) ) {
            message(paste("Cannot find array", i))
            next;
        }
        
        sectionName <- as.character(BLData@sectionData$Targets[i,"sectionName"]);
        cat(sectionName, "\n");
        
        locs <- obtainLocs(fileName = BLData@sectionData$Targets$greenLocs[i], filePath = BLData@sectionData$Targets$directory[i]);

        comb <- BeadDataPackR:::combineFiles(BLData[[i]][,c("ProbeID", "Grn", "GrnX", "GrnY")], locs)
        comb <- comb[order(comb[,5]),]
        
        for(j in 1:nSegs) {
            seg <- comb[(((j-1) * beadsPerSeg) + 1):(j * beadsPerSeg),];
            meanVar <- scoreRegistration(seg);
            res[[paste(sectionName, "Segment", j, "Grn", sep = "_")]] <- meanVar;
			## records the entire segment from the locs
			coords[[paste(sectionName, "Segment", j, "Grn", sep = "_")]] <- seg[, 3:4];
            ## records the corner coordinates of the segment
            corners[[paste(sectionName, "Segment", j, "Grn", sep = "_")]] <- seg[c(beadsPerSeg, nRows, beadsPerSeg - nRows + 1, 1), 3:4]
            ## p95 value
            p95 <- c(p95, quantile(seg[,2], probs = 0.95, names = FALSE, na.rm = TRUE))
            tiffs <- c(tiffs, file.path(BLData@sectionData$Targets[i,"directory"], BLData@sectionData$Targets[i,"greenImage"]))
			metrics <- rbind(metrics, BLData@sectionData$Metrics[i,])
        }
        
        ## detect red channel and check registration if appropriate
        if("Red" %in% colnames(BLData[[1]])) {
            
            locs <- obtainLocs(fileName = BLData@sectionData$Targets$redLocs[i], filePath = BLData@sectionData$Targets$directory[i]);
            comb <- BeadDataPackR:::combineFiles(BLData[[i]][,c("ProbeID", "Red", "RedX", "RedY")], locs)
            comb <- comb[order(comb[,5]),]
            
            for(j in 1:nSegs) {
                seg <- comb[(((j-1) * beadsPerSeg) + 1):(j * beadsPerSeg),];
                meanVar <- scoreRegistration(seg);
                res[[paste(sectionName, "Segment", j, "Red", sep = "_")]] <- meanVar;
				coords[[paste(sectionName, "Segment", j, "Red", sep = "_")]] <- seg[, 3:4];
                corners[[paste(sectionName, "Segment", j, "Red", sep = "_")]] <- seg[c(beadsPerSeg, nRows, beadsPerSeg - nRows + 1, 1), 3:4]
                p95 <- c(p95, quantile(seg[,2], probs = 0.95, names = FALSE, na.rm = TRUE))
                tiffs <- c(tiffs, file.path(BLData@sectionData$Targets[i,"directory"], BLData@sectionData$Targets[i,"redImage"]))
            }
        }
        
    }
    regScores <- new(Class = "beadRegistrationData");
    regScores@layout$twoColor <- any(grepl("Red", colnames(BLData[[1]])))
    regScores@layout$nSections <- length(array);
    regScores@layout$nSegs <- nSegs;
    regScores@registrationData <- res;
    regScores@coordinateData <- coords;
    regScores@cornerData <- corners;
    regScores@p95 <- p95;
    regScores@imageLocations <- tiffs;
	if(!is.null(metrics))
	    regScores@metrics <- metrics;
    return(regScores);
}



# checkCorners <- function(regScores, selection = NULL) {
# 
#     if(is.null(selection))
#         selection <- menu(names(regScores@registrationData), graphics = FALSE, "Selection")
# 
#     tiff <- readTIFF(regScores@imageLocations[selection]);
# 
#     idx <- strsplit(regScores@imageLocations[selection], .Platform$file.sep)[[1]]
#     idx <- idx[length(idx)]
#     col <- ifelse(grepl("Grn", idx), "green", "red")
# 
#     par(mfrow = c(2,2), mar = c(2,2,2,2))
# 
#     plotTIFF(tiff, 
#     xrange = c(max(regScores@cornerData[[selection]][1,1] - 50, 0), max(regScores@cornerData[[selection]][1,1] - 50, 0) + 300),
#     yrange = c(min(abs(regScores@cornerData[[selection]][1,2] - 250), nrow(tiff) - 301), min(abs(regScores@cornerData[[selection]][1,2] + 50), nrow(tiff) - 1)) )
#     points(regScores@coordinateData[[selection]][,1:2], pch = 19, col = col, cex = 0.3)
# 
#     plotTIFF(tiff, 
#     xrange = c(min(abs(regScores@cornerData[[selection]][2,1] + 50), ncol(tiff) - 1) - 300, min(abs(regScores@cornerData[[selection]][2,1] + 50), ncol(tiff) - 1)), 
#     yrange = c(min(abs(regScores@cornerData[[selection]][1,2] - 250), nrow(tiff) - 301), min(abs(regScores@cornerData[[selection]][1,2] + 50), nrow(tiff) - 1)) )
#     points(regScores@coordinateData[[selection]][,1:2], pch = 19, col = col, cex = 0.3)
# 
# 
#     plotTIFF(tiff, 
#     xrange = c(max(regScores@cornerData[[selection]][3,1] - 50, 0), max(regScores@cornerData[[selection]][3,1] - 50, 0) + 300), 
#     yrange = c(max(regScores@cornerData[[selection]][3,2] - 50, 0), max(regScores@cornerData[[selection]][3,2] - 50, 0) + 300) )
#     points(regScores@coordinateData[[selection]][,1:2], pch = 19, col = col, cex = 0.3)
# 
#     plotTIFF(tiff, 
#     xrange = c(min(abs(regScores@cornerData[[selection]][4,1] + 50), ncol(tiff) - 1) - 300, min(abs(regScores@cornerData[[selection]][4,1] + 50), ncol(tiff) - 1)), 
#     yrange = c(max(regScores@cornerData[[selection]][4,2] - 50, 0), max(regScores@cornerData[[selection]][4,2] - 50, 0) + 300) )
#     points(regScores@coordinateData[[selection]][,1:2], pch = 19, col = col, cex = 0.3)
# 
# }

            
            

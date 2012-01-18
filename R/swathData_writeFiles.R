writeOutFiles<-function(Swaths, an="arrayname", textstring=".txt", fullOutput = TRUE, twocolour = FALSE, outputDir = NULL){
    
    ##################################################
    ## Function to write new "bead-level" text files
    ##
    ## Arguments: Swaths - List containing two, 6-column, matrices (one for each swath)
    ##                  Columns: 1-4 Standard bead-level, 5 is an index of overlapped beads
    ##                  bead pairs in the two images have the same value 
    
    ## Output:  FULL - beads in both images are writen to a bead-level text file, with a 5th column
    ##          indicating a weight.  Beads that appear in both swaths have weight 0.5, otherwise weight 1
    ##          BASIC - overlapping beads are removed from Swath2 and no weights are output
    ###################################################

    S1 <- Swaths[[1]]
    S2 <- Swaths[[2]]
    if(length(Swaths) == 3)
        S3 <- Swaths[[3]]
    
    S1 <- S1[order(S1[,1],S1[,ncol(S1)]),]
    S2 <- S2[order(S2[,1],S2[,ncol(S2)]),]
    if(length(Swaths) == 3)
        S3 <- S3[order(S3[,1],S3[,ncol(S3)]),]
        
    ## remove the non-decoded beads
    S1 <- S1[-which(S1[,1] == 0),]
    S2 <- S2[-which(S2[,1] == 0),]
    if(length(Swaths) == 3)
        S3 <- S3[-which(S3[,1] == 0),]

    ## round coordinates so they match Illumina's format
    S1[,3:4] <- .Call("roundLocsFileValues", S1[,3:4], PACKAGE = "BeadDataPackR");
    S2[,3:4] <- .Call("roundLocsFileValues", S2[,3:4], PACKAGE = "BeadDataPackR");
    if(length(Swaths) == 3)
        S3[,3:4] <- .Call("roundLocsFileValues", S3[,3:4], PACKAGE = "BeadDataPackR");
    if(twocolour) {
        S1[,6:7] <- .Call("roundLocsFileValues", S1[,6:7], PACKAGE = "BeadDataPackR");
        S2[,6:7] <- .Call("roundLocsFileValues", S2[,6:7], PACKAGE = "BeadDataPackR");
        if(length(Swaths) == 3)
            S3[,6:7] <- .Call("roundLocsFileValues", S3[,6:7], PACKAGE = "BeadDataPackR");
    }

    if(fullOutput){
        S1[S1[,ncol(S1)] == 1, ncol(S1)] <- 0.5
        S1[S1[,ncol(S1)] == 0, ncol(S1)] <- 1
        colnames(S1)[ncol(S1)] <- "Weight"
        
        S2[S2[,ncol(S2)] == 1, ncol(S2)] <- 0.5
        S2[S2[,ncol(S2)] == 0, ncol(S2)] <- 1
        colnames(S2)[ncol(S2)] <- "Weight"
        
        if(length(Swaths) == 3) {
            S3[S3[,ncol(S3)] == 1, ncol(S3)] <- 0.5
            S3[S3[,ncol(S3)] == 0, ncol(S3)] <- 1
            colnames(S3)[ncol(S3)] <- "Weight"
        }
    }
    else {
        ## if we don't want full output remove the weights
        ## and the overlapping beads from the second swath
        S1 <- S1[ ,1:(ncol(S1) - 1)]
        S2 <- S2[which(S2[,ncol(S2)] == 0), 1:(ncol(S2) - 1)]
        if(length(Swaths) == 3)
            S3 <- S3[which(S3[,ncol(S3)] == 0), 1:(ncol(S3) - 1)]
    }
    
    write.table(S1, file=file.path(outputDir, paste(an, "-Swath1",textstring, sep = "")), row.names=FALSE, sep="\t", quote = FALSE)
    write.table(S2, file=file.path(outputDir, paste(an, "-Swath2",textstring, sep = "")), row.names=FALSE, sep="\t", quote = FALSE)
    if(length(Swaths) == 3)
        write.table(S3, file=file.path(outputDir, paste(an, "-Swath3",textstring, sep = "")), row.names=FALSE, sep="\t", quote = FALSE)
}


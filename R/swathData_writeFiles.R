writeOutFiles<-function(Swaths, an="arrayname", textstring=".txt", fullOutput = TRUE, twocolour = FALSE){

    S1<-Swaths[[1]]
    S2<-Swaths[[2]]

    ## round coordinates so they match Illumina's format
    S1[,3:4] <- .Call("roundLocsFileValues", S1[,3:4], PACKAGE = "BeadDataPackR");
    S2[,3:4] <- .Call("roundLocsFileValues", S2[,3:4], PACKAGE = "BeadDataPackR");
    if(twocolour) {
        S1[,6:7] <- .Call("roundLocsFileValues", S1[,6:7], PACKAGE = "BeadDataPackR");
        S2[,6:7] <- .Call("roundLocsFileValues", S2[,6:7], PACKAGE = "BeadDataPackR");
    }

    if(fullOutput){

        S1<-S1[order(S1[,1],S1[,ncol(S1)]),]
        S2<-S2[order(S2[,1],S2[,ncol(S2)]),]

        S1[S1[,ncol(S1)] > 0, ncol(S1)]<-0.5
        S1[S1[,ncol(S1)] == 0, ncol(S1)]<-1
        S2[S2[,ncol(S2)] > 0, ncol(S2)]<-0.5
        S2[S2[,ncol(S2)] == 0, ncol(S2)]<-1

        S1<-S1[,-(ncol(S1)-1)]
        S2<-S2[,-(ncol(S2)-1)]

        write.table(S1,row.names=F,sep="\t",file=paste(an, "-Swath1",textstring, sep = ""), quote = FALSE)
        write.table(S2,row.names=F,sep="\t",file=paste(an, "-Swath2",textstring, sep = ""), quote = FALSE)

    }

    else {

        S1<-S1[order(S1[,1]),1:(dim(S1)[2]-1)]

        S2<-S2[which(S2[,dim(S2)[2]]==0),]
        S2<-S2[order(S2[,1]),1:(dim(S2)[2]-1)]

        write.table(S1,row.names=F,sep="\t",file=paste(an, "-Swath1",textstring, sep = ""), quote = FALSE)
        write.table(S2,row.names=F,sep="\t",file=paste(an, "-Swath2",textstring, sep = ""), quote = FALSE)

    }

}


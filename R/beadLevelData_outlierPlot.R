
outlierplot <- function(BLData, array = 1, transFun = logGreenChannelTransform, outlierFun = illuminaOutlierMethod, n=3,wtsname=NULL,horizontal = TRUE, nSegments = NULL, lowOutlierCol = "blue", highOutlierCol = "pink", outlierPch = ".", main="",...){
    
    ##Find all outliers on the array
        
    wts<-1
    if(!is.null(wtsname)){wts<-getBeadData(BLData,array=array,what=wtsname)}
    
    oList = outlierFun(transFun(BLData, array=array), probeList = getBeadData(BLData,array=array,what="ProbeID"),wts = wts,n=n,...)

    cat(length(oList), " outliers found on the section\n")

    beadMeans = lapply(split(getBeadData(BLData, what="Grn", array=array), getBeadData(BLData, what="ProbeID", array=array)), mean,na.rm=TRUE)

    resids = getBeadData(BLData, what="Grn",array=array) - unlist(beadMeans[match(getBeadData(BLData, what="ProbeID", array=array), names(beadMeans))])

    oCols = NULL

    oCols[which(resids > 0 )] = highOutlierCol
    oCols[which(resids < 0 )] = lowOutlierCol
    oCols[which(resids == 0 )] = lowOutlierCol


    plotBeadLocations(BLData, array=array, BeadIDs = oList, horizontal = horizontal, col=oCols,pch=outlierPch,main=main)

    ##Work out the position of segments

    ## can we read the sdf file?
    sdfFileName <- file.path(BLData@sectionData$Targets$directory[1], list.files(as.character(BLData@sectionData$Targets$directory[1]), pattern = ".sdf")[1]);
    if( file.exists(sdfFileName) && (is.null(nSegments)) ){
        sdf <- simpleXMLparse(readLines(sdfFileName, warn = FALSE))
        nSegments <- as.integer(sdf$RegistrationParameters$SizeBlockY[[1]]);
    }
    
    ## only plot the segments if we have a value and this is a BeadChip rather than a SAM
    if( !is.null(nSegments) && !is.null(BLData@experimentData$platformClass) && !grepl("Matrix", BLData@experimentData$platformClass) ) {
        ys = getBeadData(BLData, what="GrnY", array=array)
        ys = ys - min(ys)
                
        segEnds = seq(from=0, to = max(ys), by = max(ys)/(nSegments))

        if(horizontal) 
            abline(v=segEnds, lty=2, col="red")
        else 
            abline(h=segEnds, lty=2, col="red") 
    }

}



	
calculateOutlierStats = function(BLData, array=array, transFun = logGreenChannelTransform, outlierFun = illuminaOutlierMethod, n=3, useLocs = TRUE, nSegments = 9,... ){

##Find all outliers on the array
oList = outlierFun(transFun(BLData, array=array), probeList = BLData[[array]][,1],n=n,...)

cat(length(oList), " outliers found on the section\n")

beadMeans = lapply(split(getBeadData(BLData, what="Grn", array=array), getBeadData(BLData, what="ProbeID", array=array)), mean,na.rm=TRUE)

resids = getBeadData(BLData, what="Grn",array=array) - unlist(beadMeans[match(getBeadData(BLData, what="ProbeID", array=array), names(beadMeans))])


ys = getBeadData(BLData, what="GrnY", array=array)

ys = ys - min(ys)

segEnds = seq(from=0, to = max(ys), by = max(ys)/(nSegments))

##find out which segment each bead belongs to

segInfo = cut(ys, segEnds)

table(segInfo[oList]) / table(segInfo)*100

}




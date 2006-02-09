"createBeadSummaryData" <-
function(BLData, log=FALSE, n=3, arrays=seq(1:length(BLData$R[1,]))){

#Check the object is of class BeadLevelList
  if(class(BLData) != "BeadLevelList"){
    stop("BeadLevelList object required!")
  }

#library(limma)

  probes = sort(unique(BLData$ProbeID[BLData$ProbeID[,1] > 0,1]))

  len = length(arrays)

  R = G = Rb = Gb = beadstdev = nobeads = ProbeID = nooutliers = matrix(nrow = length(probes), ncol=len)
 
  for(i in 1:length(arrays)){

    intProbeID <- as.integer(BLData$ProbeID[,arrays[i]])

    print(i)

    for(j in 1:length(probes)){
      o = findBeadStatus(BLData, probes[j], arrays[i], log=log, n=n, outputValid = TRUE, intProbeID=intProbeID)

      R[j,i] = round(mean(BLData$R[o$valid,arrays[i]],na.rm=TRUE),3)

      if(!is.null(BLData$G)){
        G[j,i] = round(mean(BLData$G[o$valid,arrays[i]],na.rm=TRUE),3)
      }
   
Rb[j,i] = round(mean(BLData$Rb[o$valid,arrays[i]],na.rm=TRUE),3)

   if(!is.null(BLData$Gb)){
     Gb[j,i] = round(mean(BLData$Gb[o$valid,arrays[i]],na.rm=TRUE),3)
   }
   
beadstdev[j,i]  = round(sd(BLData$R[o$valid,arrays[i]], na.rm=TRUE),3)

nobeads[j,i] = length(o$valid)

ProbeID[j,i] = probes[j]

nooutliers[j,i] = length(o$outliers)

}

}

BSData = list()

BSData$R = R
BSData$G = G
BSData$Rb = Rb
BSData$Gb = Gb
BSData$beadstdev = beadstdev
BSData$nobeads = nobeads
BSData$ProbeID = ProbeID
BSData$nooutliers = nooutliers
BSData$SAMPLE = BLData$SAMPLE

BSData = new("BeadSummaryList", BSData)
}


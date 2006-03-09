"createBeadSummaryData" <-
function(BLData, log=FALSE, n=3, ignoreList=NULL, arrays=seq(1:length(BLData$R[1,]))){

#Check the object is of class BeadLevelList
  if(class(BLData) != "BeadLevelList"){
    stop("BeadLevelList object required!")
  }

#library(limma)

  noprobes = length(unique(BLData$ProbeID[,1]))


  len = length(arrays)

  R = G = Rb = Gb = beadstdev = nobeads = nooutliers = matrix(0,nrow = noprobes, ncol=len)
 
  for(i in 1:length(arrays)){

    probes=sort(unique(BLData$ProbeID[,arrays[i]]))

    probes=probes[probes>0 & !is.na(probes)]
    
   
    intProbeID <- as.integer(BLData$ProbeID[,arrays[i]])

    print(i)

    for(j in 1:length(probes)){
      print(j)
      o = findBeadStatus(BLData, probes[j], arrays[i], log=log, n=n, outputValid = TRUE, intProbeID=intProbeID, ignoreList=ignoreList)

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

nooutliers[j,i] = length(o$outliers)

}

}

BSData = list()

BSData$R = R
BSData$G = G
BSData$Rb = Rb
BSData$Gb = Gb
BSData$BeadStDev = beadstdev
BSData$Nobeads = nobeads
BSData$ProbeID = unique(BLData$ProbeID[,1])
BSData$nooutliers = nooutliers
BSData$SAMPLE = BLData$SAMPLE

class(BSData) = "BeadSummaryList"

test = list()
  
if(!is.null(design)){

  ids = unique(design)

  for(i in 1:length(ids)){

    reps= which(design == ids[i])

    test[[i]] = rbind(BSData[,reps[1]], BSData[,reps[2]])
 }

}


  
}

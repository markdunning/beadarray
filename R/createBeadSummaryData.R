"createBeadSummaryData" <-
function(BLData, log=FALSE, n=3, ignoreList=NULL, arrays=seq(1:length(BLData$R[1,])), design = rep(1, ncol(BLData))){

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
     
   
    intProbeID <- as.integer(sort(BLData$ProbeID[,arrays[i]]))
    probeIndex <- c(1:length(intProbeID))
    probeIndex <- probeIndex[sort.list(BLData$ProbeID[,arrays[i]])]
    
    print(i)
    nextStart = 1
#    o = sapply(probes, findBeadStatus, BLData = BLData, array = arrays[i], log=log, n=n, outputValid = TRUE, intProbeID=intProbeID, ignoreList=ignoreList, probeIndex = probeIndex)
    
    for(j in 1:length(probes)){

#      o = findBeadStatus(BLData, probes[j], array = arrays[i], log=log, n=n, outputValid = TRUE, intProbeID=intProbeID[c(nextStart:length(intProbeID))], ignoreList=ignoreList, probeIndex = probeIndex[c(nextStart:length(intProbeID))])

      o = findBeadStatus(BLData, probes[j], array = arrays[i], log=log, n=n, outputValid = TRUE, intProbeID=intProbeID, ignoreList=ignoreList, probeIndex = probeIndex, startSearch = nextStart)
      
      R[j,i] = mean(BLData$R[o$valid,arrays[i]],na.rm=TRUE)

      if(!is.null(BLData$G)){
        G[j,i] = mean(BLData$G[o$valid,arrays[i]],na.rm=TRUE)
      }
   
      Rb[j,i] = mean(BLData$Rb[o$valid,arrays[i]],na.rm=TRUE)

      if(!is.null(BLData$Gb)){
        Gb[j,i] = mean(BLData$Gb[o$valid,arrays[i]],na.rm=TRUE)
      }
      
      beadstdev[j,i]  = sd(BLData$R[o$valid,arrays[i]], na.rm=TRUE)

      nobeads[j,i] = length(o$valid)

      nooutliers[j,i] = length(o$outliers)
      nextStart = o$nextStart
 #     print(nextStart)
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
#  BSData$SAMPLE = BLData$SAMPLE

  class(BSData) = "BeadSummaryList"

  test = list()
  
  if(!is.null(design)){

    ids = unique(design)

    for(i in 1:length(ids)){
      cat("i: ",i,"\n")
      reps = which(design == ids[i])
      test[[i]] <- BSData[,reps[1]]

      if(length(reps) > 1){
      
        for(j in 2:length(reps)){
          test[[i]] = rbind(test[[i]], BSData[,reps[j]])
        }
      }
    }
    BSData2 <- test[[1]]
    if(length(test) > 1){
    for(k in 2:length(test)){
        BSData2 <- cbind(BSData2, test[[k]])
      }
    }
    p1 = sort(unique(BLData$ProbeID[,1]))
    p2 = sort(unique(BLData$ProbeID[,2]))
    
    BSData2$ProbeID = union(p1,p2)
    BSData2
  }

  else{

    BSData

  }
  
}

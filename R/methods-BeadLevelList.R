setMethod("initialize", "BeadLevelList",
          function(.Object,
									
									 G = new("matrix"),
									 Gb = new("matrix"),
									 GrnX = new("matrix"),
									 GrnY= new("matrix"),
									 ProbeID = new("matrix"),
									 targets = new("data.frame")
									) {
									.Object@G<-G
									.Object@Gb<-Gb
									.Object@GrnX<-GrnX
									.Object@GrnY<-GrnY
									.Object@ProbeID<-ProbeID
									.Object@targets<-targets
									
									
									
									
								  .Object
})





findBeadStatus <- function(BLData, probes, array, log=FALSE, n=3,
                           outputValid = FALSE, intProbeID = NULL, ignoreList=NULL,
                           probeIndex = NULL, startSearch = 1){

  if(is.null(intProbeID)){
    intProbeID <- as.integer(sort(BLData$ProbeID[,array]))
    probeIndex <- c(1:length(intProbeID))
    probeIndex <- probeIndex[sort.list(BLData$ProbeID[,array])]
  }
 
  outliers = valid = vector()
  
  for(i in 1:length(probes)){
    if(! (probes[i] %in% ignoreList)){
  
      temp = getProbeIndicesC(BLData, probe = probes[i], intProbe = intProbeID,
        index = probeIndex, startSearch = startSearch)
      probe_ids = temp[[1]]
      startSearch = temp[[2]]

      inten <- BLData$G[probe_ids,array]
  
   #nas will be a list of beads which have NA intensity
    nas=NULL

    if(length(which(is.na(inten)))>0){

      nas = probe_ids[is.na(inten) | inten < 0]
      probe_ids = probe_ids[!is.na(inten)]
      inten = inten[!is.na(inten)]
    }

    if(log){
      raw_inten = log2(BLData$G[probe_ids,array])
      raw_inten = raw_inten[!is.na(raw_inten)]
    }
    else{
      raw_inten = inten
    }
  
    m = mean(inten, na.rm=TRUE, trim=0.5)
    ma = mad(inten, na.rm=TRUE)
  index = (inten > m + n *ma | inten < m - n*ma | raw_inten < 0)
  
    outliers.temp=probe_ids[index]

    outliers =c(outliers, outliers.temp, nas)
    if(outputValid)
      valid = c(valid, probe_ids[!index])
  }
  }
  if(outputValid){
    result <- list(outliers = outliers, valid = valid, nextStart = startSearch)
    result
  }
  else
    outliers
}

getProbeIndicesC <- function(BLData, probe, intProbe, index, startSearch = 1){

  if(is.null(BLData$ProbeID)) stop("ProbeID column was not found in BLData object")

  ind <- .C("findIndices", as.integer(probe), intProbe, as.integer(nrow(BLData)), result = integer(length = 25000),
            pos = as.integer(startSearch), PACKAGE="beadarray")

  
  ind2 <- vector()
  i = 1;
  while(ind$result[i] != 0){
    ind2  <- c(ind2, index[ind$result[i]])
    i = i+1;
  }
  list(ind2, ind$pos)
}




findAllOutliers <- function(BLData, array, log = FALSE, n = 3){

  probes <- sort(unique(BLData$ProbeID[BLData$ProbeID[,array] > 0,array]))

  if(log){
    finten <- log2(BLData$G[,array])
    finten[is.na(finten)] = 0
  }
  else{
    finten <- BLData$G[,array]
  }
  probeList <- BLData$ProbeID[,array]
  nbeads <- length(BLData$G[,array])

  start = 0

  foo <- .C("findAllOutliers", as.double(finten), binStatus = integer(length = nbeads), as.integer(probeList), as.integer(probes), as.integer(length(probes)), as.integer(nbeads), as.integer(start), PACKAGE = "beadarray")

  which((probeList > 0) & (foo$binStatus == 0))

}

"getProbeIntensities" <-
function(BLData, ProbeIDs, array,log=TRUE){

#Check the object is of class BeadLevelList
  if(class(BLData) != "BeadLevelList"){
    stop("BeadLevelList object required!")
  }
  
if(log){
log2(BLData$G[BLData$ProbeID[,array] %in% ProbeIDs,array])
}
else{
BLData$G[BLData$ProbeID[,array] %in% ProbeIDs,array]
}

}




createBeadSummaryData <- function(BLData, log = FALSE, n = 3, arrays=seq(1:length(BLData$G[1,])), imagesPerArray = 2, probes = NULL){

  #Check the object is of class BeadLevelList
  if(class(BLData) != "BeadLevelList"){
    stop("BeadLevelList object required!")
  }

  len = ncol(BLData)

  if(imagesPerArray == 1){
    temp <- BLData[BLData$ProbeID[,1] != 0,1]
  }
  else if(imagesPerArray == 2){
    temp <- rbind(BLData[BLData$ProbeID[,1] != 0,1], BLData[BLData$ProbeID[,2] != 0,2])
  }
  else{
    stop("You can only specify 1 or 2 images per array")
  }
  
  if(is.null(probes)){
    probes = sort(unique(as.vector(temp$ProbeID)))
  }
    probes = probes[probes>0 & !is.na(probes)]
    noprobes = length(probes)

  if(imagesPerArray == 1){
    R = G = Rb = Gb = BeadStDev = NoBeads = nooutliers = matrix(0,nrow = noprobes, ncol=len) }
  else{
     R = G = Rb = Gb = BeadStDev = NoBeads = nooutliers = matrix(0,nrow = noprobes, ncol=(len/2)) }

  i = j = 1
   while(j <= len){
    print(i)
    if(log){
     finten <- log2(temp$G)
     binten <- log2(temp$Gb)
    }
    else {
      finten <- temp$G
      binten <- temp$Gb
    }
     probeIDs <- as.integer(temp$ProbeID)
#    start = (length(which($ProbeID[,i] == 0)))x
     start = 0
     blah <- .C("createBeadSummary",  as.double(finten),  as.double(binten), probeIDs, as.integer(probes), as.integer(noprobes), as.integer(length(temp$G)),
                 fore = double(length = noprobes), back = double(length = noprobes), sd = double(length = noprobes), noBeads = integer(length = noprobes),
                 noOutliers = integer(length = noprobes), nextStart = as.integer(start), PACKAGE = "beadarray")

     G[,i] = blah$fore
     Gb[,i] = blah$back
     BeadStDev[,i] = blah$sd
     NoBeads[,i] = blah$noBeads
     nooutliers[,i] = blah$noOutliers
     j = j+imagesPerArray
     i = i + 1
     rm(probeIDs, blah)
     gc()
     if((imagesPerArray == 1) && (i <= ncol(BLData))){
       temp = BLData[BLData$ProbeID[,i] != 0, i]
     }
     else if((imagesPerArray == 2) && (j < ncol(BLData))){
       temp = rbind(BLData[BLData$ProbeID[,j] != 0,j], BLData[BLData$ProbeID[,j+1] != 0,j+1])
     }


   }

##Change the standard deviation to the standard error

BeadStDev = BeadStDev / sqrt(NoBeads)

BSData<-new("ExpressionSetIllumina")

assayData(BSData)=assayDataNew(exprs = G, BeadStDev=BeadStDev, NoBeads = NoBeads)
rownames(exprs(BSData)) = probes
 
BSData
}






readBeadSummaryData <- function(targets=NULL, header=T, sep=",",path=NULL,
                                   columns = list(probeID = "ProbeID",
                                     AvgSig = "AvgSig", Nobeads = "Nobeads",
                                     BeadStDev = "BeadStDev", Sample = "SAMPLE",
                                     Detection = "Detection", SAM = "SAM", Avebkg="Avebkg"),
                                   other.columns = NULL)
{


if(!is.null(targets)) {
targets = as.character(targets[,1])

if(!is.null(path)) targets = file.path(path, targets)

}
else{

if(!is.null(path)){

	if(is.null(targets)) {targets = setwd(path);targets=dir()
	}

	else{
	targets=dir()
	}

}

 }

  filecounter = arraycounter = 0
  r = read.table(targets[1], header = header, sep = sep)
  cat("Reading file ", targets[1], "\n")
  cat("First line", "\n")
  print(r[1,])

  #calculate how many probes are on each array
  nrows = length(unique(r[,columns$probeID]))
  #calculate how many arrays are in these files
  ArraysPerFile = ((length(r[,columns$probeID]))/nrows)
  if(ArraysPerFile == 1){
    R = G = Rb = Gb = beadstdev = nobeads = probeID = nooutliers = Detection =
      matrix(nrow = nrows, ncol=(length(targets)*ArraysPerFile))
    sample = SAM = vector(mode = "character", length = ArraysPerFile*(length(targets)))

    #	Other columns
	if(!is.null(other.columns)) {
		other.columns <- as.character(other.columns)
		j <- match(other.columns,colnames(r),0)
		if(any(j>0)) {
			other.columns <- colnames(r)[j]
			other <- list()
			for (j in other.columns) other[[j]] <- matrix(NA, nrow = nrow(R), ncol = ncol(R)) 
		} else {
			other.columns <- NULL
		}
	}
    
    for(i in 1:length(targets)){
      if(i == 1){
        r = r
      }
      else{
        r = read.table(targets[i], header = header, sep = sep)
        cat("Reading file ", targets[i], "\n")
  	  cat("First line", "\n")
  	  print(r[1,])
}
      for(j in 1:ArraysPerFile){
#        d = r[(((j-1)*nrows)+1):(j*nrows),]
#        cat("J:",j,"\n")
#        cat("Count:",filecounter,"\n")
        probeID[,j+filecounter] = r[,columns$probeID]
        R[,j+filecounter] = r[,columns$AvgSig]
        Rb[,j+filecounter] = r[,columns$Avebkg]
        beadstdev[,j+filecounter] = r[,columns$BeadStDev]
        nobeads[,j+filecounter] = r[,columns$Nobeads]
        Detection[,j+filecounter] = r[,columns$Detection]
        sample[j+filecounter] = as.character(r[1,columns$Sample])
        SAM[j+filecounter] = as.character(r[1,columns$SAM])
        if(!is.null(other.columns))
          for (m in other.columns) {
            other[[m]][,j+filecounter] <- r[,j] 
          }
       }
       filecounter = filecounter + ArraysPerFile
    }
    BSData = list()
    BSData$R = R
    BSData$Rb = Rb
    BSData$probeID = probeID
    BSData$beadstdev = beadstdev
    BSData$nobeads = nobeads
    BSData$Detection = Detection
    BSData$SAMPLE = sample
    BSData$SAM = SAM
    if(!is.null(other.columns))
      BSData$other = other
    BSData = new("BeadSummaryList", BSData)
  }
  else {
    for(i in 1:length(targets)){
      if(i == 1){
        r = r
        R = G = Rb = Gb = beadstdev = nobeads = probeID = nooutliers = Detection =
          matrix(nrow = nrows, ncol=ArraysPerFile)
        sample = SAM = vector(mode = "character", length = ArraysPerFile)
#reading in any values for "other.columns"        
        if(!is.null(other.columns)) {
          other.columns <- as.character(other.columns)
          j <- match(other.columns,colnames(r),0)
          if(any(j>0)) {
            other.columns <- colnames(r)[j]
            other <- list()
            for (j in other.columns) other[[j]] <- matrix(NA, nrow = nrow(R), ncol = ncol(R)) 
          } else {
            other.columns <- NULL
          }
	}
        
        for(j in 1:ArraysPerFile){
          d = r[(((j-1)*nrows)+1):(j*nrows),]
          
          probeID[,j+filecounter] = d[,columns$probeID]
          R[,j+filecounter] = d[,columns$AvgSig]
          Rb[,j+filecounter] = d[,columns$Avebkg]
          beadstdev[,j+filecounter] = d[,columns$BeadStDev]
          nobeads[,j+filecounter] = d[,columns$Nobeads]
          Detection[,j+filecounter] = d[,columns$Detection]
          sample[j+filecounter] = as.character(d[1,columns$Sample])
          SAM[j+filecounter] = as.character(d[1,columns$SAM])
          if(!is.null(other.columns))
          for (m in other.columns) {
            other[[m]][,j+filecounter] <- d[,m] 
          }
        } 
        filecounter = filecounter + ArraysPerFile
      }
      
      else{
        r = read.table(targets[i], header = header, sep = sep)
        cat("Reading file ", targets[1], "\n")
  	  cat("First line", "\n")
        print(r[1,])
 
        ArraysPerFile = ((length(r[,columns$probeID]))/nrows)
         R.temp = G.temp = Rb.temp = Gb.temp = beadstdev.temp = nobeads.temp = probeID.temp = nooutliers.temp = Detection.temp =
          matrix(nrow = nrows, ncol=ArraysPerFile)
        sample.temp = SAM.temp = vector(mode = "character", length = ArraysPerFile)
#reading in any values for "other.columns"
        if(!is.null(other.columns)) {
          other.columns <- as.character(other.columns)
          j <- match(other.columns,colnames(r),0)
          if(any(j>0)) {
            other.columns <- colnames(r)[j]
            other.temp <- list()
            for (j in other.columns) other.temp[[j]] <- matrix(NA, nrow = nrow(R.temp), ncol = ncol(R.temp)) 
          } else {
            other.columns <- NULL
          }
	}
        
      for(j in 1:ArraysPerFile){
        d = r[(((j-1)*nrows)+1):(j*nrows),]
        
        probeID.temp[,j] = d[,columns$probeID]
        R.temp[,j] = d[,columns$AvgSig]
        Rb.temp[,j] = d[,columns$Avebkg]
        beadstdev.temp[,j] = d[,columns$BeadStDev]
        nobeads.temp[,j] = d[,columns$Nobeads]
        Detection.temp[,j] = d[,columns$Detection]
        sample.temp[j] = as.character(d[1,columns$Sample])
        SAM.temp[j] = as.character(d[1,columns$SAM])
        if(!is.null(other.columns))
          for (m in other.columns) {
            other.temp[[m]][,j] <- d[,m] 
          }
      }
        filecounter = filecounter + ArraysPerFile
        probeID <- cbind(probeID, probeID.temp)
        R <- cbind(R,R.temp)
        Rb<-cbind(Rb, Rb.temp)
        beadstdev <- cbind(beadstdev, beadstdev.temp)
        nobeads <- cbind(nobeads, nobeads.temp)
        Detection <- cbind(Detection, Detection.temp)
        sample <- append(sample, sample.temp)
        SAM <- append(SAM, SAM.temp)
        if(!is.null(other.columns)){
          for(k in 1:length(other)){
            other[[k]] <- cbind(other[[k]],other.temp[[k]])
          }
        }
      }
    BSData = list()
    BSData$R = R
    BSData$Rb = Rb
    BSData$probeID = probeID
    BSData$beadstdev = beadstdev
    BSData$nobeads = nobeads
    BSData$Detection = Detection
    BSData$SAMPLE = sample
    BSData$SAM = SAM
    if(!is.null(other.columns))
      BSData$other = other
    BSData = new("BeadSummaryList", BSData)
    }
  }
}

#readBeadSummaryData <- function(files){
#  r = read.table(files[1], header=T, sep=",")
#  R = G = Rb = Gb = beadstdev = nobeads = probeID = nooutliers = matrix(nrow = nrow(r), ncol=length(files))
#  sample = vector(length = length(files))
#  for(i in 1:length(files)){
#    r = read.table(files[i], header=T, sep=",")
#    print(dim(r))
#    probeID[,i] = r[,1]
#    R[,i] = r[,2]
#    beadstdev[,i] = r[,3]
#    nobeads[,i] = r[,4]
#    sample[i] = r[1,6]
#  }
#  BSData = list()
#  BSData$R = R
#  BSData$probeID = probeID
#  BSData$beadstdev = beadstdev
#  BSData$nobeads = nobeads
#  BSData = new("BeadSummaryList", BSData)
#}


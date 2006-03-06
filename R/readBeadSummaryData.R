readBeadSummaryData <- function(targets=NULL, header=T, sep=",",path=NULL,
                                   columns = list(ProbeID = "ProbeID",
                                     AvgSig = "AvgSig", Nobeads = "Nobeads",
                                     BeadStDev = "BeadStDev", Detection="Detection"),
                                   other.columns = NULL)
{


if(!is.null(targets)) {
targets = as.character(targets[,1])

if(!is.null(path)) targets = file.path(path, targets)

}


else{

if(!is.null(path)){

setwd(path)
}

targets=dir()


}




  filecounter = arraycounter = 0
  r = read.table(targets[1], header = header, sep = sep)
  cat("Reading file ", targets[1], "\n")
  cat("First line", "\n")
  print(r[1,])

  #calculate how many probes are on each array
  nrows = length(unique(r[,columns$ProbeID]))
  
   #calculate how many arrays are in these files
  
#Find how many column names match the names we expect

avg.cols = which(match(strtrim(colnames(r), nchar(columns$AvgSig)), columns$AvgSig)==1)
nobeads.cols = which(match(strtrim(colnames(r), nchar(columns$Nobeads)), columns$Nobeads)==1)
beadstdev.cols = which(match(strtrim(colnames(r), nchar(columns$BeadStDev)), columns$BeadStDev)==1)
detection.cols = which(match(strtrim(colnames(r), nchar(columns$Detection)), columns$Detection)==1)


#How many columns match the column heading for the Signal. ie how many arrays are there
matches = length(avg.cols)


#By default we expect each array to be listed in columns going across the file
readAcross = TRUE
#If only one column match was found then array are listed under each other in the file


  #We assume that arrays are   
  if (matches >1) {
    ArraysPerFile = matches
  }

  else{
    readAcross=FALSE
    ArraysPerFile = ((length(r[,columns$ProbeID]))/nrows)
  
    avg.cols = rep(avg.cols, ArraysPerFile)
    nobeads.cols = rep(nobeads.cols, ArraysPerFile)
    beadstdev.cols = rep(beadstdev.cols, ArraysPerFile)
    detection.cols = rep(detection.cols, ArraysPerFile)

  }


    

    for(i in 1:length(targets)){


      if(i==1){

        
        r = read.table(targets[i], header = header, sep = sep)
        cat("Reading file ", targets[i], "\n")
  	cat("First line", "\n")
  	print(r[1,])
        R =  beadstdev = nobeads = ProbeID = nooutliers = Detection =
        matrix(nrow = nrows, ncol=ArraysPerFile)


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

      for(j in 1:ArraysPerFile){

        if(readAcross){
          rowsToRead = 1:nrows
        }
        else{
          rowsToRead = (((j-1)*nrows)+1):(j*nrows)
        }
        ProbeID[,j+filecounter] = as.character(r[rowsToRead,columns$ProbeID])
        R[,j+filecounter] = r[rowsToRead,avg.cols[j]]
        beadstdev[,j+filecounter] = r[rowsToRead,beadstdev.cols[j]]
        nobeads[,j+filecounter] = r[rowsToRead,nobeads.cols[j]]
        Detection[,j+filecounter] = r[rowsToRead,detection.cols[j]]
        
        if(!is.null(other.columns)){
          for (m in other.columns) {
            other[[m]][,j+filecounter] <- r[,j] 
          }
       }
       
      }

      filecounter = filecounter + ArraysPerFile  
      }

      
    else{
      r = read.table(targets[i], header = header, sep = sep)
       cat("Reading file ", targets[i], "\n")
  	cat("First line", "\n")
        print(r[1,])

       if(!readAcross){
      
        ArraysPerFile = ((length(r[,columns$ProbeID]))/nrows)
        avg.cols = rep(which(match(strtrim(colnames(r), nchar(columns$AvgSig)), columns$AvgSig)==1),ArraysPerFile)
         nobeads.cols = rep(which(match(strtrim(colnames(r), nchar(columns$Nobeads)), columns$Nobeads)==1),ArraysPerFile)
         beadstdev.cols = rep(which(match(strtrim(colnames(r), nchar(columns$BeadStDev)), columns$BeadStDev)==1),ArraysPerFile)
         detection.cols = rep(which(match(strtrim(colnames(r), nchar(columns$Detection)), columns$Detection)==1),ArraysPerFile)
      }

       else{
         avg.cols = which(match(strtrim(colnames(r), nchar(columns$AvgSig)), columns$AvgSig)==1)
         nobeads.cols = which(match(strtrim(colnames(r), nchar(columns$Nobeads)), columns$Nobeads)==1)
         beadstdev.cols = which(match(strtrim(colnames(r), nchar(columns$BeadStDev)), columns$BeadStDev)==1)
         detection.cols = which(match(strtrim(colnames(r), nchar(columns$Detection)), columns$Detection)==1)

         

         ArraysPerFile = length(avg.cols) 
       }
         
      
         R.temp  = beadstdev.temp = nobeads.temp = ProbeID.temp =  Detection.temp =
          matrix(nrow = nrows, ncol=ArraysPerFile)

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
        if(readAcross){
          rowsToRead = 1:nrows
        }
        else{
          rowsToRead = (((j-1)*nrows)+1):(j*nrows)
        }
        
        
        
        ProbeID.temp[,j] =as.character(r[rowsToRead,columns$ProbeID])
        R.temp[,j] = r[rowsToRead,avg.cols[j]]
        beadstdev.temp[,j] = r[rowsToRead,beadstdev.cols[j]]
        nobeads.temp[,j] =r[rowsToRead,nobeads.cols[j]]
        Detection.temp[,j] = r[rowsToRead,detection.cols[j]]


        if(!is.null(other.columns))
          for (m in other.columns) {
            other.temp[[m]][rowsToRead,j] <- r[,m] 
          }
        
        

      }


    
      ProbeID <- cbind(ProbeID, ProbeID.temp)
        R <- cbind(R,R.temp)
      
        beadstdev <- cbind(beadstdev, beadstdev.temp)
        nobeads <- cbind(nobeads, nobeads.temp)
        Detection <- cbind(Detection, Detection.temp)

        if(!is.null(other.columns)){
          for(k in 1:length(other)){
            other[[k]] <- cbind(other[[k]],other.temp[[k]])
          }
        }
    }


    filecounter = filecounter + ArraysPerFile
   
    }
    BSData = list()
    BSData$R = R
   
    BSData$ProbeID = ProbeID
    BSData$beadstdev = beadstdev
    BSData$nobeads = nobeads
    BSData$Detection = Detection
    if(!is.null(other.columns)){
      BSData$other = other
    }
    BSData = new("BeadSummaryList", BSData)
    BSData$targets = targets
    BSData
}




readBeadSummaryData<- function(targets=NULL, header=T, sep=",",path=NULL,
                                   columns = list(ProbeID = "TargetID",
                                     AvgSig = "AVG_Signal", Nobeads = "Avg_NBEADS",
                                     BeadStDev = "BEAD_STDEV", Detection="Detection"),
                                   other.columns = NULL, skip = 7)
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

  r = try(read.table(targets[1], header = header, sep = sep, skip = skip, nrows = 1))
  cat("Reading file ", targets[1], "\n")

  temp <- scan(file = targets[1], skip = skip+1, quiet=TRUE,sep = sep, what = as.list(rep("character",  ncol(r))))


  #cat("First line", "\n")
  #print(r[1,])

  #calculate how many probes are on each array
  nrows = length(temp[[1]])
  
   #calculate how many arrays are in these files
  
#Find how many column names match the names we expect

   
#How many columns match the column heading for the Signal. ie how many arrays are there

#By default we expect each array to be listed in columns going across the file
readAcross = TRUE
#If only one column match was found then array are listed under each other in the file

avg.cols = which(match(strtrim(colnames(r), nchar(columns$AvgSig)), columns$AvgSig)==1)

matches = length(avg.cols)

if(length(which(match(strtrim(colnames(r), nchar(columns$ProbeID)), columns$ProbeID)==1))==0) stop(paste("Could not find any columns called", columns$ProbeID))



  #We assume that arrays are   
  if (matches >1) {
    ArraysPerFile = matches
  }

  else{
    readAcross=FALSE
    ArraysPerFile = nrow(r) / length(unique(r[,columns$ProbeID]))
  }


    BSData = list()

    for(i in 1:length(targets)){
      r = read.table(targets[i], header = header, sep = sep, skip = skip)
       cat("Reading file ", targets[i], "\n")   

      if(i==1){
        if(readAcross){
#First item in columns is the ProbeID columnd, we only read this once 
        BSData = readCols(r, columns=columns[2:5])
        

    #	Other columns
	if(!is.null(other.columns)) {

        BSData$other = readCols(r, columns=other.columns, colnames=other.columns) 
	}

          
      }

   
     else{
       for(j in 1:ArraysPerFile){
         if(j==1){

           
         
         r= r[(((j-1)*nrows)+1):(j*nrows),]

         BSData = readCols(r, columns=columns[2:5])
         
         if(!is.null(other.columns)) {

           BSData$other = readCols(r, columns=other.columns, colnames=other.columns) 
         }

       }

         else{
          r= r[(((j-1)*nrows)+1):(j*nrows),]

         BSData = cbind(BSData, readCols(r, columns=columns[2:5]))
         
         if(!is.null(other.columns)) {

           BSData$other = cbind(BSData$other,readCols(r, columns=other.columns, colnames=other.columns)) 
         }
           
       }
         
     }
     }
    }   
    else{
      
      if(readAcross){
      
      temp = readCols(r, columns=columns[2:5])
         
      
        if(!is.null(other.columns)) {
          temp$other = readCols(r, columns=other.columns, colnames=other.columns)
	}


    }
      else{

        
       for(j in 1:ArraysPerFile){

         
         r= r[(((j-1)*nrows)+1):(j*nrows),]

         temp= readCols(r, columns=columns[2:5])
         if(!is.null(other.columns)) {

           temp$other = readCols(r, columns=other.columns, colnames=other.columns) 
         }
       }
         
       }

      
    BSData = cbind(BSData, temp)
      
        
    }

    }
    BSData$ProbeID = as.character(unique(r[,columns$ProbeID]))
    BSData = new("BeadSummaryList", BSData)
    BSData$targets = targets
    BSData
}


readCols <- function(table, columns = list(AvgSig = "AVG_Signal", Nobeads = "Avg_NBEADS",
                              BeadStDev = "BEAD_STDEV", Detection="Detection"),
                     colnames = c("R", "Nobeads", "BeadStDev", "Detection")){
  
  cols.to.read=list(length=length(columns))

  nrows=nrow(table)
  
for(i in 1:length(columns)){

cols.to.read[[i]] = which(match(strtrim(colnames(table), nchar(columns[i])), columns[i])==1)

if(length(cols.to.read[[i]])==0) stop(paste("Could not find any columns called", columns[i]))

}



  #How many columns match the column heading for the Signal. ie how many arrays are there


  output=list()
  
  for(i in 1:length(columns)){

    output[[i]] = table[, cols.to.read[[i]]]

  }

  names(output) <- colnames
  output

  output = new("BeadSummaryList", output)

}



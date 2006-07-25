"readBeadLevelData" <-
  function(targets, path=NULL, imageManipulation = "none", backgroundSize=17, columns=list(ProbeID="Code", x="x", y="y"), numrow = NULL){
	
   if(is.null(targets$Image1)) stop("Error: beadTargets object must contain Image1 column")
   if(is.null(targets$xyInfo)) stop("Error: beadTargets object must contain xyInfo column")

   manip = switch(imageManipulation, none = 0, sharpen = 1, sasik = 2, 3)

   if(manip == 3){
     stop("The imageManipulation arguement must be one of: \"none\", \"sharpen\" or \"sasik\"")
   }
   
   #Take this line out later and make it optional
   foregroundCalc = "Illumina"

   if(foregroundCalc == "Illumina"){
     fground = 0
   }
   else if(foregroundCalc == "sasik"){
     fground = 1
   }
     
	
   tifFiles = as.character(targets$Image1)
   csv_files = as.character(targets$xyInfo)

   if(!is.null(targets$Image2)) tifFiles2 = as.character(targets$Image2)   			
  	
   if(is.null(path)) path=getwd()	

    tifs = dir(path=path, pattern =".tif")
    csvs = dir(path=path, pattern =".csv")

   #needs to be changed for two color data
    for(i in 1:length(csv_files)){
      if(! (csv_files[i] %in% csvs)){
        stop(paste("File not found ", csv_files[i]))
      }
      if(! (tifFiles[i] %in% tifs)){
        stop(paste("File not found ", tifFiles[i]))
      }  
    }
 
    k = nrow(targets)

   if(is.null(numrow)){
     csvNcol = ncol(read.table(csv_files[1], sep=",", header=T, nrows = 1))
     if(csvNcol == 4){ #One colour data
       r = scan(csv_files[[1]], what = list(integer(0), NULL, NULL, NULL), sep = ",", skip = 1, quiet = TRUE)
     }
     else if(csvNcol == 6){ #Two colour data
       r = scan(csv_files[[1]], what = list(integer(0), NULL, NULL, NULL, NULL, NULL), sep = ",", skip = 1, quiet = TRUE)
     }
     else { #Stop if there is a weird number of columns i.e. not 4 or 6
       stop("Unknown format for xyFile")
     }
     numrow = length(r[[1]])
   }
     
   BLData <- list(R = matrix(nrow = numrow, ncol=k), Rb = matrix(nrow = numrow, ncol=k),
                  x = matrix(nrow = numrow, ncol=k), y = matrix(nrow = numrow, ncol=k),
                  ProbeID = matrix(nrow = numrow, ncol=k))
  
    if(!is.null(targets$Image2)){
      G = Gb = matrix(nrow = numrow, ncol=k)
    }

   ##requires large changes to read 2 colour data
   for(i in 1:k){

     file=csv_files[i]
     
     if(!is.null(path)) file=file.path(path, file) 
     
     dat1 <- scan(file, what = list(ProbeID = integer(0), NULL, x = numeric(0), y = numeric(0)), sep = ",", skip = 1, quiet = TRUE)
     ord <- order(dat1$ProbeID)
     
     BLData$x[,i] <- dat1$x[ord]
     BLData$y[,i] <- dat1$y[ord]
     BLData$ProbeID[,i] <- dat1$ProbeID[ord]
     
     rm(dat1)
     
     numBeads = length(BLData$x[,i])
     
     intensities <- .C("readBeadImage", as.character(tifFiles[i]), as.double(BLData$x[,i]),
                       as.double(BLData$y[,i]), as.integer(numBeads), foreGround = double(length = numBeads),
                       backGround = double(length = numBeads), as.integer(backgroundSize), as.integer(manip),
                       as.integer(fground), PACKAGE = "beadarray")

     BLData$R[,i] <- intensities[[5]]
     BLData$Rb[,i] <- intensities[[6]]

#     if(csvNcol == 6){
#       redIntensities <- .C("readBeadImage", as.character(tifFiles[i]), as.double(BLData$x[,i]),
#                            as.double(BLData$y[,i]), as.integer(numBeads), foreGround = double(length = numBeads),
#                            backGround = double(length = numBeads), as.integer(backgroundSize), as.integer(manip),
#                            as.integer(fground), PACKAGE = "beadarray")
#     }
     
   }

   BLData$targets = targets
   
   BLData$backgroundSize = backgroundSize

   BLData$normalised = 0
   
   BLData$backgroundCorrected = 0
   
   class(BLData) = "BeadLevelList"
   BLData
   
 }                                                              

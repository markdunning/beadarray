"readBeadLevelData" <-
  function(targets, path=NULL, imageManipulation = "sharpen", backgroundSize=17, columns=list(ProbeID="Code", x="x", y="y"), sep = "\t", numrow = NULL){
	
   if(is.null(targets$Image1)) stop("Error: beadTargets object must contain Image1 column")
   if(is.null(targets$xyInfo)) stop("Error: beadTargets object must contain xyInfo column")

   manip = switch(imageManipulation, none = 0, sharpen = 1, sasik = 2, sasikFaster = 3, 4)

   if(manip == 4){
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
  	
   if(is.null(path)) path=getwd()	

#    tifs = dir(path=path, pattern =".tif")
#    csvs = dir(path=path, pattern =".csv")

   #Check that the specified xyFiles and images exist
   files = dir(path = path)
   for(i in 1:length(csv_files)){
     if(!(csv_files[i] %in% files)){
       stop(paste("File not found ", csv_files[i]))
     }
     if(! (tifFiles[i] %in% files)){
        stop(paste("File not found ", tifFiles[i]))
      }
   }

   ##Check if there is a second image specified check they exist
   if(!is.null(targets$Image2)){
     tifFiles2 = as.character(targets$Image2)
     for(i in 1:length(csv_files)){
       if(!(tifFiles2[i] %in% files)){
         stop(paste("File not found ", tifFiles2[i]))
       }
     }
   }
   
   #needs to be changed for two color data
#    for(i in 1:length(csv_files)){
#      if(! (csv_files[i] %in% csvs)){
#        stop(paste("File not found ", csv_files[i]))
#      }
#      if(! (tifFiles[i] %in% tifs)){
#        stop(paste("File not found ", tifFiles[i]))
#      }  
#    }
 
    k = nrow(targets)

   if(is.null(numrow)){
     csvNcol = ncol(read.table(csv_files[1], sep=sep, header=T, nrows = 1))
     if(csvNcol == 4){ #One colour data
       r = scan(csv_files[[1]], what = list(integer(0), NULL, NULL, NULL), sep = sep, skip = 1, quiet = TRUE)
     }
     else if(csvNcol == 7){ #Two colour data
       r = scan(csv_files[[1]], what = list(integer(0), NULL, NULL, NULL, NULL, NULL, NULL), sep = sep, skip = 1, quiet = TRUE)
     }
     else { #Stop if there is a weird number of columns i.e. not 4 or 6
       stop("Unknown format for xyFile")
     }
     numrow = length(r[[1]])
   }

   if(csvNcol == 4){
     BLData <- list(G = matrix(nrow = numrow, ncol=k), Gb = matrix(nrow = numrow, ncol=k),
                    GrnX = matrix(nrow = numrow, ncol=k), GrnY = matrix(nrow = numrow, ncol=k),
                    ProbeID = matrix(nrow = numrow, ncol=k))
   }
   else if(csvNcol == 7){
     BLData <- list(G = matrix(nrow = numrow, ncol=k), Gb = matrix(nrow = numrow, ncol=k),
                    R = matrix(nrow = numrow, ncol=k), Rb = matrix(nrow = numrow, ncol=k),
                    GrnX = matrix(nrow = numrow, ncol=k), GrnY = matrix(nrow = numrow, ncol=k),
                    RedX = matrix(nrow = numrow, ncol=k), RedY = matrix(nrow = numrow, ncol=k),
                    ProbeID = matrix(nrow = numrow, ncol=k))
   }
  
#    if(!is.null(targets$Image2)){
#      R = Rb = matrix(nrow = numrow, ncol=k)
#    }

   for(i in 1:k){

     file=csv_files[i]
     
     if(!is.null(path)) file=file.path(path, file) 

     if(csvNcol == 4){
       dat1 <- scan(file, what = list(ProbeID = integer(0), NULL, GrnX = numeric(0),
                            GrnY = numeric(0)), sep = sep, skip = 1, quiet = TRUE)
     }
     else if(csvNcol == 7){
       dat1 <- scan(file, what = list(ProbeID = integer(0), NULL, GrnX = numeric(0),
                            GrnY = numeric(0), NULL, RedX = numeric(0), RedY = numeric(0)),
                    sep = sep, skip = 1, quiet = TRUE)
     }
     ord <- order(dat1$ProbeID)
     
     BLData$GrnX[,i] <- dat1$GrnX[ord]
     BLData$GrnY[,i] <- dat1$GrnY[ord]
     BLData$ProbeID[,i] <- dat1$ProbeID[ord]

     if(csvNcol == 7){
       BLData$RedX[,i] <- dat1$RedX[ord]
       BLData$RedY[,i] <- dat1$RedY[ord]
     }
     
     rm(dat1)
     gc()
     
     numBeads = length(BLData$GrnX[,i])
     
     greenIntensities <- .C("readBeadImage", as.character(tifFiles[i]), as.double(BLData$GrnX[,i]),
                       as.double(BLData$GrnY[,i]), as.integer(numBeads), foreGround = double(length = numBeads),
                       backGround = double(length = numBeads), as.integer(backgroundSize), as.integer(manip),
                       as.integer(fground), PACKAGE = "beadarray")

     BLData$G[,i] <- greenIntensities[[5]]
     BLData$Gb[,i] <- greenIntensities[[6]]

     rm(greenIntensities)
     gc()

     if(csvNcol == 7){
       redIntensities <- .C("readBeadImage", as.character(tifFiles2[i]), as.double(BLData$RedX[,i]),
                            as.double(BLData$RedY[,i]), as.integer(numBeads), foreGround = double(length = numBeads),
                            backGround = double(length = numBeads), as.integer(backgroundSize), as.integer(manip),
                            as.integer(fground), PACKAGE = "beadarray")

       BLData$R[,i] <- redIntensities[[5]]
       BLData$Rb[,i] <- redIntensities[[6]]
     }
     
   }

   BLData$targets = targets
   
   BLData$backgroundSize = backgroundSize

   BLData$normalised = 0
   
   BLData$backgroundCorrected = 0
   
   class(BLData) = "BeadLevelList"
   BLData
   
 }                                                              

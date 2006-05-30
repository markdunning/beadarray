"readBeadImagesSlow" <-
  function(targets, path=NULL, sharpen=TRUE, backgroundSize=17, columns=list(ProbeID="Code", x="x", y="y"), numrow = NULL){
	
   if(is.null(targets$Image1)) stop("Error: beadTargets object must contain Image1 column")
   if(is.null(targets$xyInfo)) stop("Error: beadTargets object must contain xyInfo column")
	
   pgm_files1 = as.character(targets$Image1)
   csv_files = as.character(targets$xyInfo)
   if(!is.null(targets$Image2)) pgm_files2 = as.character(targets$Image2)
   			
  	
    if(is.null(path)) path=getwd()	

    pgms = dir(path=path, pattern =".pgm")
    csvs = dir(path=path, pattern =".csv")

    for(i in 1:length(csv_files)){

      if(! (csv_files[i] %in% csvs)){

        stop(paste("File not found ", csv_files[i]))

      }

      if(! (pgm_files1[i] %in% pgms)){

        stop(paste("File not found ", pgm_files1[i]))

      }

      
    }
 
    k = nrow(targets)

   if(is.null(numrow)){
#     r=read.table(csv_files[1], sep=",", header=T)
     r = scan(csv_files[[1]], what = list(integer(0), NULL, NULL, NULL, NULL, NULL), sep = ",", skip = 1, quiet = TRUE)
     numrow = length(r[[1]])
   }
     
    #Read a set of k arrays using images and csv files
    
 #   R = Rb =  x = y = ProbeID = matrix(nrow = numrow, ncol=k)

   BLData <- list(R = matrix(nrow = numrow, ncol=k), Rb = matrix(nrow = numrow, ncol=k),
                  x = matrix(nrow = numrow, ncol=k), y = matrix(nrow = numrow, ncol=k),
                  ProbeID = matrix(nrow = numrow, ncol=k))

  
    if(!is.null(targets$Image2)){

      G = Gb = matrix(nrow = numrow, ncol=k)

    }

    for(i in 1:k){

      file=csv_files[i]

	if(!is.null(path)) file=file.path(path, file) 

      dat1 <- scan(file, what = list(NULL, ProbeID = integer(0),NULL, NULL, x = numeric(0), y = numeric(0)), sep = ",", skip = 1, quiet = TRUE)
      
      BLData$x[,i] <- dat1$x
      BLData$y[,i] <- dat1$y
      BLData$ProbeID[,i] <- dat1$ProbeID
#      dat = read.table(file, sep=",", header=T)  
#      xs = dat[,columns$x] + 1
#      ys = dat[,columns$y] + 1
#      BLData$ProbeID[,i] = dat[,columns$ProbeID]
#      BLData$x[,i] = xs
#      BLData$y[,i] = ys

      rm(dat1)
      cat("Calculating foreground intensities for", pgm_files1[i], "\n")

       file=pgm_files1[i]

	if(!is.null(path)) file=file.path(path, file)       

 	I = read.pgmfile(file, sep=",", header=T)

#Must add 1 to the x and y coordinates

      BLData$R[,i] = calculateForegroundIntensities(I, BLData$x[,i]+1, BLData$y[,i]+1, sharpen=sharpen)
      cat("Calculating background intensities.\n")
      BLData$Rb[,i] = calculateBackground(I, BLData$x[,i]+1, BLData$y[,i]+1, n=backgroundSize)
      
      if(!is.null(targets$Image2)){
	  file=pgm_files2[i]

	  if(!is.null(path)) file=file.path(path, file)   	

        I = read.pgmfile(file, sep=",", header=T)
        
        G[,i] = calculateForegroundIntensities(I, BLData$xs[,i], BLData$ys[,i], sharpen=sharpen)

        Gb[,i] = calculateBackground(I, BLData$xs[,i], BLData$ys[,i], n=backgroundSize)              
        }
      rm(I)
#      rm(xs)
#      rm(ys)
#      gc()
    }
#    BLData = list()
#    BLData$R = R

#BLData$Rb = Rb

if(!is.null(targets$Image2)){

BLData$G = G
BLData$Gb = Gb

}

#   BLData$x = x

#   BLData$y = y

#   BLData$ProbeID = ProbeID

    BLData$targets = targets

    BLData$sharpened = as.numeric(sharpen)

    BLData$backgroundSize = backgroundSize

    BLData$normalised = 0

    BLData$backgroundCorrected = 0

#    BLData = new("BeadLevelList", BLData)
   class(BLData) = "BeadLevelList"
   BLData
  }


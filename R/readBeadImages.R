"readBeadImages" <-
  function(targets, path=NULL, sharpen=TRUE, backgroundSize=17, storeRawData=FALSE, columns=list(probeID="Code", x="x", y="y")){
	
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
    r=read.table(csv_files[1], sep=",", header=T)
     
    #Read a set of k arrays using images and csv files
    
    R = Rb =  x = y = probeID = matrix(nrow = nrow(r), ncol=k)

  
    if(!is.null(targets$Image2)){

      G = Gb = matrix(nrow = nrow(r), ncol=k)

    }

    if(storeRawData){

      rawR = rawG = matrix(nrow=nrow(r), ncol=k)

    }

	rm(r)
gc()

    for(i in 1:k){


      file=csv_files[i]

	if(!is.null(path)) file=file.path(path, file) 


      dat = read.table(file, sep=",", header=T)
    
      xs = dat[,columns$x] + 1
      ys = dat[,columns$y] + 1


      probeID[,i] = dat[,columns$probeID]

      x[,i] = xs

      y[,i] = ys

      rm(dat)
      cat("Calculating foreground intensities for", pgm_files1[i], "\n")

       file=pgm_files1[i]

	if(!is.null(path)) file=file.path(path, file)       

 	I = read.pgmfile(file, sep=",", header=T)

      R[,i] = calculateForegroundIntensities(I, xs, ys, sharpen=sharpen)
      cat("Calculating background intensities.\n")
      Rb[,i] = calculateBackground(I, xs, ys, n=backgroundSize)

      if(storeRawData){

        rawR[,i] = calculateForegroundIntensities(I, xs, ys, sharpen=FALSE)

      }

      if(!is.null(targets$Image2)){
	  file=pgm_files2[i]

	  if(!is.null(path)) file=file.path(path, file)   	

        I = read.pgmfile(file, sep=",", header=T)


        
        G[,i] = calculateForegroundIntensities(I, xs, ys, sharpen=sharpen)

        Gb[,i] = calculateBackground(I, xs, ys, n=backgroundSize)
        
        if(storeRawData){
  
          rawG[,i] = calculateForegroundIntensities(I, xs, ys, sharpen=FALSE)
          
        }
        rm(I)
      }
    }
    BLData = list()
    BLData$R = R

BLData$Rb = Rb

if(!is.null(targets$Image2)){

BLData$G = G

BLData$Gb = Gb

}

BLData$x = x

BLData$y = y

BLData$probeID = probeID

if(storeRawData){

BLData$other$rawR = rawR

if(!is.null(targets$Image2)){

BLData$other$rawG = rawG

}

}
    BLData$targets = targets

    BLData$sharpened = as.numeric(sharpen)

    BLData$backgroundSize = backgroundSize

    BLData$normalised = 0

    BLData$backgroundCorrected = 0

    BLData = new("BeadLevelList", BLData)


  }


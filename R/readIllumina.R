"readIllumina" =
  function(arrayNames=NULL, path=".", textType=".txt", 
           annoPkg=NULL, beadInfo=NULL, useImages=TRUE, 
           singleChannel=TRUE, targets=NULL, 
           imageManipulation = "sharpen", backgroundSize=17,
           storeXY=TRUE, sepchar="_", dec=".", metrics=FALSE,
           metricsFile="Metrics.txt", backgroundMethod="subtract",
           offset=0, normalizeMethod="none", tiffExtGrn="_Grn.tif", 
           tiffExtRed="_Red.tif",...){

  if(textType==".csv") sep=","
  else sep="\t"	

  if(!useImages) {
    message("'useImages=FALSE': please check that the 'singleChannel'
            argument is appropriate for your data set", sep="")
    if(is.null(arrayNames))
       message(paste("It is advisable to specify 'arrayNames' to avoid reading in spurious", textType, "files from the specified directory"))
  }
  xyFiles = dir(path=path, pattern =textType)
  metricpos = grep(metricsFile, xyFiles)
  if(length(metricpos)>0) {
     xyFiles = xyFiles[-metricpos]
  }
  if(length(xyFiles)==0)
    stop(paste("No files with extension ", textType, "found"))
  arrays = strtrim(xyFiles, nchar(xyFiles)-4)  
  if(!is.null(arrayNames)) {
    # check whether file extensions are present - if they are remove them
    arrayNames = gsub(textType, "", arrayNames)
    arrays = arrayNames[which(arrayNames %in% arrays)] 
    if(length(arrays)==0)      
      stop("'arrayNames' did not match with files in this directory")
    }
  if(useImages) {
    manip = switch(imageManipulation, none = 0, sharpen = 1, morph = 2, morphFaster = 3, 4)

    if(manip == 4){
      stop("The imageManipulation argument must be one of: \"none\", \"sharpen\" or \"morph\"")
    }
   
    # Take this line out later and make it optional
    foregroundCalc = "Illumina"

    if(foregroundCalc == "Illumina"){
      fground = 0
    }
    else if(foregroundCalc == "morph"){
      fground = 1
    }
     
    GImages = dir(path=path, pattern =tiffExtGrn)
    RImages = dir(path=path, pattern =tiffExtRed)
                                        
    if(length(GImages)==0)
      stop("No tiffs found")

    ###Find which files have both Green Images and xy information
    arrays = intersect(arrays, strtrim(GImages, nchar(GImages)-8))	
  
    ##Check to see if we have two channels
    if(length(RImages)!=0) {
      arrays = intersect(arrays, strtrim(RImages, nchar(RImages)-8))
      if(singleChannel==TRUE)
        cat("Red images found and 'singleChannel=TRUE'.  Only Green images will be read in.\n")
    }
    else { # no Red images found
      if(singleChannel==FALSE) {
        cat("No red images found, setting singleChannel=TRUE\n")
        singleChannel = TRUE
      }
    }
  }
  
  k = length(arrays)
  cat("Found", k, "arrays","\n")
   
  BLData = new("BeadLevelList")


####Check the specified annotation name

if(is.null(annoPkg)) warning("No annotation package was specified. Need to use SetAnnotation later\n")

else{
     data(ExpressionControlData)
     m=match(annoPkg, names(ExpressionControlData))
    
     if(is.na(m)) warning("Annotation must be one of ", paste(names(ExpressionControlData)," "))
     else BLData@annotation = annoPkg	
}




  usedIDs=NULL
  arrayInfo=list(arrayNames=as.character(arrays), 
             nBeads=rep(0, length(arrays)), 
             background=backgroundMethod,
             normalization=normalizeMethod)

  if(length(grep(sepchar, arrays))==length(arrays)) {
     tmp = unlist(strsplit(arrays, sepchar))
     if(length(tmp)%%3==0) {
        arrayInfo$chip = tmp[seq(1,length(tmp), by=3)]
        arrayInfo$row = tmp[seq(2,length(tmp), by=3)]
        arrayInfo$col = tmp[seq(3,length(tmp), by=3)]
     }
     else {
        arrayInfo$chip = tmp[seq(1,length(tmp), by=2)]
        arrayInfo$row = tmp[seq(2,length(tmp), by=2)]
        arrayInfo$col = rep(1, length(arrayInfo$chip))
     }
   }                                                                                                       

  csv_files = vector(length=length(arrays))
  csv_files = file.path(path, paste(arrays, textType,sep=""))
  csvNcol = ncol(read.table(csv_files[1], sep=sep, dec=dec, header=T, nrows = 1))

  if(singleChannel) {
    arrayInfo$channels = "single"
    if(storeXY){
      ncolumns = 6
      headings=c("ProbeID", "G", "Gb","GrnX", "GrnY", "wts")	
    }
    else{
      ncolumns = 4
      headings=c("ProbeID", "G", "Gb", "wts")
    }
  }
  else { # two-channel data
    arrayInfo$channels = "two"
    tifFiles2 = vector(length=length(arrays))
    tifFiles2 = file.path(path, paste(arrays, tiffExtRed, sep=""))
    if(storeXY){
      ncolumns = 8
      headings=c("ProbeID", "G", "Gb", "GrnX", "GrnY", "R", "Rb", "wts")
    }
    else{
      ncolumns = 6
      headings=c("ProbeID", "G", "Gb","R", "Rb", "wts")
      }
  }

  if(useImages) {
    tifFiles = vector(length=length(arrays))
    tifFiles = file.path(path, paste(arrays, tiffExtGrn, sep=""))

    for(i in 1:k) {
   	 
     file=csv_files[i]
     
     if(csvNcol == 4){
       fc = file(file, open="r")
       dat1 = scan(file=fc, what = list(ProbeID = integer(0), 
                   NULL, GrnX = numeric(0), GrnY = numeric(0)), 
                   sep = sep, dec=dec, skip = 1, quiet = TRUE)
       close(fc)
     }


     else if(csvNcol == 7){
       fc = file(file, open="r")
       dat1 = scan(file=fc, what = list(ProbeID = integer(0), 
                   NULL, GrnX = numeric(0), GrnY = numeric(0), 
                   NULL, RedX = numeric(0), RedY = numeric(0)),
                   sep = sep, dec=dec, skip = 1, quiet = TRUE)
       close(fc)                   
     }

     ##An older single channel format 
     else if(csvNcol==6){
       fc = file(file, open="r")
       dat1 = scan(file=fc, what = list(NULL,
                   ProbeID = integer(0),NULL, NULL, 
                   GrnX = numeric(0), GrnY = numeric(0)), 
                   sep = sep, dec=dec, skip = 1, quiet = TRUE)
       close(fc) 
       RedX = dat1$GrnX
       RedY = dat1$GrnY
     }
     
     ord = order(dat1$ProbeID)
     data = matrix(nrow=length(ord), ncol = ncolumns)
     colnames(data) = headings	
     
     data[,1] = dat1$ProbeID[ord]
     usedIDs = union(usedIDs, unique(data[,1]))

     numBeads = length(dat1$GrnX)
     arrayInfo$nBeads[i] = numBeads
     
     greenIntensities = .C("readBeadImage", as.character(tifFiles[i]), as.double(dat1$GrnX[ord]),
                       as.double(dat1$GrnY[ord]), as.integer(numBeads), foreGround = double(length = numBeads),
                       backGround = double(length = numBeads), as.integer(backgroundSize), as.integer(manip),
                       as.integer(fground),PACKAGE = "beadarray")

     l = length(unique(data[,1]))
#     e = .C("startEndPos", as.integer(data[,1]),as.integer(numBeads),integer(length=as.integer(l)), integer(length=as.integer(l)))
#     endPos[[i]]=e[[4]]
     
     cat("Background correcting: method =", backgroundMethod, "\n")
     data[,2] = bgCorrectSingleArray(fg=greenIntensities[[5]], bg=greenIntensities[[6]], method=backgroundMethod, offset=offset)
     data[,3] = greenIntensities[[6]]

	if(storeXY){
     data[,4] = (dat1$GrnX[ord] - min(dat1$GrnX))
     data[,5] = (dat1$GrnY[ord] - min(dat1$GrnY))
	}

     rm(greenIntensities)
     gc()

     if(!singleChannel){
       if(csvNcol ==7) RedX=dat1$RedX;RedY=dat1$RedY
       
       redIntensities = .C("readBeadImage", as.character(tifFiles2[i]), as.double(RedX[ord]),
                            as.double(RedY[ord]), as.integer(numBeads), foreGround = double(length = numBeads),
                            backGround = double(length = numBeads), as.integer(backgroundSize), as.integer(manip),
                            as.integer(fground), PACKAGE = "beadarray")
       rm(dat1)
       gc()             
       if(storeXY) {
       cat("Background correcting: method =", backgroundMethod, "\n")
       data[,6] = bgCorrectSingleArray(fg=redIntensities[[5]], bg=redIntensities[[6]], method=backgroundMethod, offset=offset)
       data[,7] = redIntensities[[6]]
       }
       else {
       cat("Background correcting: method =", backgroundMethod, "\n")
       data[,4] = bgCorrectSingleArray(fg=redIntensities[[5]], bg=redIntensities[[6]], method=backgroundMethod, offset=offset) #redIntensities[[5]]
       data[,5] = redIntensities[[6]]
       }
       
       rm(redIntensities)
       gc()
       cat("Normalizing R and G intensities: method =", normalizeMethod, "\n")
       if(storeXY)
         data[,c(2,6)] = normalizeSingleArray(data[,c(2,6)], method=normalizeMethod)
       else
         data[,c(2,4)] = normalizeSingleArray(data[,c(2,4)], method=normalizeMethod)
     }
     data[,which(colnames(data) == "wts")] = 1 ##initialise wts to 1
     assign(arrays[i], as.data.frame(data), envir=BLData@beadData)
   }
 }
 else { # use intensities stored in the text files
    for(i in 1:k) {
   	 
     file=csv_files[i]
     cat("Reading raw data from", file, "\n")
     if(csvNcol == 4){
       fc = file(file, open="r")
       dat1 = scan(file=fc, what = list(ProbeID = integer(0), 
              Grn = numeric(0), GrnX = numeric(0),
              GrnY = numeric(0)), sep = sep, dec=dec, skip = 1, quiet = TRUE)
       close(fc)
     }

     else if(csvNcol == 7){
       fc = file(file, open="r")
       dat1 = scan(file=fc, what = list(ProbeID = integer(0), 
                   Grn = numeric(0), GrnX = numeric(0),
                   GrnY = numeric(0), Red = numeric(0),
                   RedX = numeric(0), RedY = numeric(0)),
                   sep = sep, dec=dec, skip = 1, quiet = TRUE)
       close(fc)                   
     }

     ##An older single channel format 
     else if(csvNcol==6){
       fc = file(file, open="r")
       dat1 = scan(file=fc, what = list(NULL,
                   ProbeID = integer(0), Grn = numeric(0), 
                   Red = numeric(0), GrnX = numeric(0), 
                   GrnY = numeric(0)), sep = sep, 
                   dec=dec, skip = 1, quiet = TRUE)
       close(fc) 
     }
     
     ord = order(dat1$ProbeID)
     data = matrix(nrow=length(ord), ncol = ncolumns)
     colnames(data) = headings	
     data[,1] = dat1$ProbeID[ord]
     usedIDs = union(usedIDs, unique(data[,1]))

     numBeads = length(dat1$GrnX)
     arrayInfo$nBeads[i] = numBeads

     data[,2]=dat1$Grn[ord]
     data[,3]=0
     cat("Background correcting: method =", backgroundMethod, "\n")
     data[,2] = bgCorrectSingleArray(fg=data[,2], bg=data[,3], method=backgroundMethod, offset=offset)

     if(storeXY) {
       if(!singleChannel) {
         data[,6]=dat1$Red[ord]
         data[,7]=0
         data[,6] = bgCorrectSingleArray(fg=data[,6], bg=data[,7], method=backgroundMethod, offset=offset)
         cat("Normalizing R and G intensities: method =", normalizeMethod, "\n")
         data[,c(2,6)] = normalizeSingleArray(data[,c(2,6)], method=normalizeMethod)
       }
       data[,4]=dat1$GrnX[ord]-min(dat1$GrnX)
       data[,5]=dat1$GrnY[ord]-min(dat1$GrnY)
     }
     else {
       if(!singleChannel) {
         data[,4]=dat1$Red[ord]
         data[,5] = 0
         data[,4] = bgCorrectSingleArray(fg=data[,4], bg=data[,5], method=backgroundMethod, offset=offset)
         cat("Normalizing R and G intensities: method =", normalizeMethod, "\n")
         data[,c(2,4)] = normalizeSingleArray(data[,c(2,4)], method=normalizeMethod)
       }
     }
     data[,which(colnames(data) == "wts")] = 1 ##initialise weights to 1
     assign(arrays[i], as.data.frame(data), envir=BLData@beadData)
   }
  }
  BLData@arrayInfo = arrayInfo


## Add targets information (if available)
if(!is.null(targets))
   BLData@phenoData = new("AnnotatedDataFrame", data=targets)
else
   BLData@phenoData = new("AnnotatedDataFrame", data=data.frame(arrayName=arrays))

## Add bead information (sequence, Illumina IDs used in annoPkg etc)
if(!is.null(beadInfo) && is.data.frame(beadInfo))
   BLData@beadAnno = beadInfo

##Look for scanner metrics file
                                           

##Create template for the various QC scores

ScanMetrics = matrix(nrow=k, ncol=13)


BASHscores = matrix(nrow=k, ncol=3)
colnames(BASHscores) = c("Number of Compact Defects","Number of Diffuse Defects", "Extended Defects Score")


outlierScores = matrix(nrow=k, ncol=9)
colnames(outlierScores) = paste("Segment", 1:9)


controlScores = matrix(nrow=k, ncol=9)

rownames(BASHscores) = rownames(ScanMetrics) = rownames(outlierScores) = rownames(controlScores) = arrays



if(metrics) {                                                          
  metrics = dir(path=path, pattern=metricsFile)
  if(length(metrics)==1)
    ScanMetrics = read.table(file.path(path, metricsFile), sep="\t", dec=dec, header=T)
}

qcScores = list("ScanMetrics"=ScanMetrics, "BASHscores"=BASHscores, "OutlierDistribution"=outlierScores, "controlProbeScores" = controlScores)
BLData@qcScores = qcScores

BLData
}                                                              


bgCorrectSingleArray = function(fg, bg, method = "subtract", offset = 0, verbose=FALSE) {
    method = match.arg(method, c("none", "subtract", "half", 
        "minimum", "edwards", "normexp", "rma"))
    switch(method,
	subtract={
            fg = fg - bg
        },
	half={
            fg = pmax(fg - bg, 0.5)
        },
	minimum={
            fg = fg - bg
            j = fg < 1e-18
            if(any(j, na.rm = TRUE)) {
              m = min(fg[!j], na.rm = TRUE)
              fg[j] = m/2
            }
        },
	edwards={
#		Log-linear interpolation for dull spots as in Edwards (2003).
#		The threshold values (delta) are chosen such that the number of
#		spots with (0 < R-Rb < delta) is f=10% of the number spots
#		with (R-Rb <= 0) for each channel and array.
#		Note slight change from Edwards (2003).
            one = matrix(1,NROW(fg),1)
	    delta.vec = function(d, f=0.1) {
	      quantile(d, mean(d<1e-16,na.rm=TRUE)*(1+f), na.rm=TRUE)
	    }
	    sub = as.matrix(fg-bg)
	    delta = one %*% apply(sub, 2, delta.vec)
	    fg = ifelse(sub < delta, delta*exp(1-(bg+delta)/fg), sub)
	},
	normexp={
	    x = fg - bg
	    out = normexp.fit(x)
            gc()
#	     if(verbose) cat("G: bg.bias=",out$par[1]," bg.sd=",exp(out$par[2])," fg.mean=",exp(out$par[3]),"\n",sep="")
	    fg = normexp.signal(out$par,x)
            gc()
#	    if(verbose) cat("Corrected array",i,"\n")
        },
        rma={
	  require("affy")
	    fg = bg.adjust(fg - bg)
        })    
	if(offset) {
            fg = fg + offset
        }
        gc()
        fg
    }

normalizeSingleArray = function(x, method = "quantile") {
    method = match.arg(method, c("none", "quantile", "vsn"))
    switch(method,
        none={
          y = x
        },
	quantile={
            y = normalizeQuantiles(x)
        },
	vsn={
            require("vsn")
            y = exprs(vsn2(x)) #@hx
        })
    y
  }

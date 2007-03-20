"readIllumina" =
  function(arrayNames=NULL, path=NULL, textType=".csv", annoFile=NULL,
           targets=NULL, imageManipulation = "sharpen", backgroundSize=17,
           storeXY=TRUE, sepchar="_", metrics=FALSE,
           metricsFile="Metrics.txt", backgroundMethod="none", offset=0,
           normalizeMethod="none", ...){

  if(textType==".csv") sep=","
  else sep="\t"	

   manip = switch(imageManipulation, none = 0, sharpen = 1, sasik = 2, sasikFaster = 3, 4)

   if(manip == 4){
     stop("The imageManipulation argument must be one of: \"none\", \"sharpen\" or \"sasik\"")
   }
   
   #Take this line out later and make it optional
   foregroundCalc = "Illumina"

   if(foregroundCalc == "Illumina"){
     fground = 0
   }
   else if(foregroundCalc == "sasik"){
     fground = 1
   }
     
   TwoChannel=FALSE
                                                                                   GImages = dir(pattern ="_Grn.tif")	
   RImages = dir(pattern ="_Red.tif")
   xyFiles = dir(pattern =textType) 

   ###Find which files have both Green Images and xy information
   arrays = intersect(strtrim(GImages, nchar(GImages)-8), strtrim(xyFiles, nchar(xyFiles)-4))	



  
   ##Check to see if we have two channels
   if(length(RImages)!=0) {

    arrays = intersect(arrays, strtrim(RImages, nchar(RImages)-8))
    TwoChannel = TRUE
   }

   if(!is.null(arrayNames))
     arrays = arrayNames[which(arrayNames %in% arrays)] 
                                                                                   csv_files = tifFiles = tifFiles2 = vector(length=length(arrays))


   csv_files = paste(arrays, textType,sep="")
   tifFiles = paste(arrays, "_Grn.tif",sep="")

  if(TwoChannel)
    tifFiles2 = paste(arrays, "_Red.tif",sep="")
   
  cat("Found", length(arrays), "arrays","\n")
  k = length(arrays)


  csvNcol = ncol(read.table(csv_files[1], sep=sep, header=T, nrows = 1))


  BLData = new("BeadLevelList")
  endPos = new("list")
  usedIDs=NULL
   if(!TwoChannel){

	if(storeXY){
     ncolumns = 5
     headings=c("ProbeID", "G", "Gb","GrnX", "GrnY")	
	}
	else{
	ncolumns=3
	headings=c("ProbeID", "G", "Gb")
	}

   }
   else {
	if(storeXY){

        ncolumns = 7
	headings=c("ProbeID", "G", "Gb", "GrnX", "GrnY", "R", "Rb")


	}
	else{
	ncolumns = 5
	headings=c("ProbeID", "G", "Gb","R", "Rb")

}

    

   }
  
#    if(!is.null(targets$Image2)){
#      R = Rb = matrix(nrow = numrow, ncol=k)
#    }

   arrayInfo=list(arrayNames=as.character(arrays), 
             nBeads=rep(0, length(arrays))) # targets=matrix(nrow=length(arrays), ncol=2)
#   targets[,1] = arrays
   
                                                                                   if(TwoChannel) arrayInfo$channels = "two"
                                                                                   else arrayInfo$channels = "single"                    

   if(length(grep(sepchar, arrays))==length(arrays)) {
     tmp = unlist(strsplit(arrays, "_"))
     arrayInfo$chip = tmp[seq(1,length(tmp), by=3)]
     arrayInfo$row = tmp[seq(2,length(tmp), by=3)]
     arrayInfo$col = tmp[seq(3,length(tmp), by=3)]
   }                                                                                                                                                         
   for(i in 1:k) {
    

     	 
     file=csv_files[i]
     
     if(!is.null(path)) file=file.path(path, file) 

     if(csvNcol == 4){
       fc = file(file, open="r")
       dat1 = scan(file=fc, what = list(ProbeID = integer(0), NULL, GrnX = numeric(0),
                            GrnY = numeric(0)), sep = sep, skip = 1, quiet = TRUE)
       close(fc)
     }


     else if(csvNcol == 7){
       fc = file(file, open="r")
       dat1 = scan(file=fc, what = list(ProbeID = integer(0), NULL, GrnX = numeric(0),
                            GrnY = numeric(0), NULL, RedX = numeric(0), RedY = numeric(0)),
                    sep = sep, skip = 1, quiet = TRUE)
       close(fc)                   
     }

     ##An older single channel format 
     else if(csvNcol==6){
       fc = file(file, open="r")
       dat1 = scan(file=fc, what = list(NULL,ProbeID = integer(0),NULL, NULL, GrnX = numeric(0),
                            GrnY = numeric(0)), sep = sep, skip = 1, quiet = TRUE)
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

     if(TwoChannel){
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
     assign(arrays[i], as.data.frame(data), envir=BLData@beadData)
   }

   BLData@arrayInfo = arrayInfo
   

#   probeindex = matrix(nrow=length(usedIDs), ncol=length(arrays))

 #  for(i in 1:length(arrays)){
     
 #    p = unique(get(arrays[i], envir=BLData@beadData)[,1])
 #    probeindex[match(p, usedIDs),i] = endPos[[i]]

 #  }
 # rownames(probeindex) = usedIDs



 # BLData@probeindex = probeindex


                                                                                ## Add bead annotation information (if available)
if(!is.null(annoFile)){
  if(!is.null(path)) annoFile=file.path(path, annoFile) 
  if(length(grep(".opa", annoFile))==1)
     BLData@beadAnno = readOPA(annoFile)
  else
     BLData@beadAnno = read.table(annoFile, ...)
}

## Add targets information (if available)
if(!is.null(targets))
   pData(BLData@phenoData) = targets

##Look for scanner metrics file
                                           
if(metrics) {                                                          
  metrics = dir(pattern=metricsFile)
  if(length(metrics)==1)
    BLData@scanMetrics = read.table(metrics, sep="\t", header=T)
}

BLData
}                                                              


#u = unique(BLData[[1]][,1])

#probeindex = matrix(nrow=length(u), ncol=length(BLData))

#for(i in 1:length(BLData)){

 # for(j in 1:1:length(u)){

  #  probeindex[j,i] = max(which(BLData[[i]][,1]==u[j]))

  #}

#}


G=function(list){

  max = 0
  p=list$probeindex

  for(i in 1:ncol(p)){
    m = max(p[,i], na.rm=TRUE)
    if(m >max) max = m
  }
  

  Gmat = matrix(nrow=max, ncol=ncol(p),0)

  for(i in 1:ncol(p)){
    Gcol = which(colnames(list[[i]]) == "G")
        if(is.null(Gcol)) stop("Could not find R column")
    Gmat[1:nrow(list[[i]]),i] = list[[i]][,Gcol]
  }

  Gmat

}

Gb=function(list){

  max = 0
  p=list$probeindex

  for(i in 1:ncol(p)){
    m = max(p[,i], na.rm=TRUE)
    if(m >max) max = m
  }
  

  Gmat = matrix(nrow=max, ncol=ncol(p),0)

  for(i in 1:ncol(p)){
    Gcol = which(colnames(list[[i]]) == "Gb")
        if(is.null(Gcol)) stop("Could not find R column")
    Gmat[1:nrow(list[[i]]),i] = list[[i]][,Gcol]
  }

  Gmat

}


R=function(list){
  max=0
  for(i in 1:length(list)){
    if(nrow(list[[i]]) > max) max = nrow(list[[i]])

  }

  Gmat = matrix(nrow=max, ncol=length(list),0)

  for(i in 1:length(list)){
    Gcol = which(colnames(list[[i]]) == "R")
        if(is.null(Gcol)) stop("Could not find R column")
    Gmat[1:nrow(list[[i]]),i] = list[[i]][,Gcol]
  }

  Gmat

}

Rb=function(list){
  max=0
  for(i in 1:length(list)){
    if(nrow(list[[i]]) > max) max = nrow(list[[i]])

  }

  Gmat = matrix(nrow=max, ncol=length(list),0)

  for(i in 1:length(list)){
    Gcol = which(colnames(list[[i]]) == "Rb")
    if(is.null(Gcol)) stop("Could not find R column")
    Gmat[1:nrow(list[[i]]),i] = list[[i]][,Gcol]
  }

  Gmat

}



Ps=function(list) {

  max = 0
  p=list$probeindex

  for(i in 1:ncol(p)) {
    m = max(p[,i], na.rm=TRUE)
    if(m >max) max = m
  }
  

  Gmat = matrix(nrow=max, ncol=ncol(p),0)

  for(i in 1:ncol(p)) {
    Gcol = which(colnames(list[[i]]) == "ProbeID")
        if(is.null(Gcol)) stop("Could not find R column")
    Gmat[1:nrow(list[[i]]),i] = list[[i]][,Gcol]
  }

  Gmat

}


readOPA <- function(OPAfile, path=NULL) { # code taken from BeadarraySNP package (methods-SnpSetIllumina.R file)
  if(!is.null(path))
    OPAfile <- file.path(path, OPAfile)
  firstfield <- scan(OPAfile, what = "", sep = ",", flush = TRUE, quiet = TRUE, blank.lines.skip = FALSE, multi.line = FALSE)
  skip <- grep("Ilmn ID", firstfield, fixed=TRUE)
  if (length(skip) == 0) stop("Cannot find \"Ilmn ID\" in OPA info file")
        enddata<- grep("[Gentrain Request]", firstfield, fixed=TRUE)
  if (length(enddata) == 0) stop("Cannot find \"[Gentrain Request]\" in OPA info file")
        OPAinfo<-read.table(OPAfile, skip=skip[1]-1, header = TRUE, sep = ",", as.is = TRUE, check.names = FALSE, nrows=enddata[1]-skip[1]-1)
  colnames(OPAinfo)<-c("Illname","snpid","oligo1","oligo2","oligo3","IllCode","IllOligo","IllStrand","snpbases","CHR","Ploidy","Species","MapInfo","TopGenomicSeq","CustomerStrand")
  rownames(OPAinfo)<-OPAinfo[,"snpid"]
  OPAinfo
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
            require("affy")
            y = normalize.quantiles(x)
        },
	vsn={
            require("vsn")
            y = vsn(x)@exprs/log(2)
        })
    y
  }

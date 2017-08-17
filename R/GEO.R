###Functions for preparing files for a GEO submission compatible with the current guidelines at;
###http://www.ncbi.nlm.nih.gov/geo/info/geo_illu.html


makeGEOSubmissionFiles <- function(normData, rawData,forceDetection=TRUE,basename="GEO",softTemplate=NULL){
  message("Creating GEO meta data...")
  createGEOMeta(normData,basename)
  message("Creating GEO processed Matrix....")
  createGEOMatrix(normData,forceDetection,basename,softTemplate,normalised=T)
  message("Creating GEO non-normalised Matrix....")
  createGEOMatrix(rawData,forceDetection,basename,softTemplate,normalised=F)
}




createGEOMeta <- function(normData,basename){
  
  data(metaTemplate)
  
  pd <- pData(normData)
  metanames <- metaTemplate[["SAMPLES"]]
  
  metaMatrix <- matrix(nrow = nrow(pd),ncol=length(metanames))
  colnames(metaMatrix) <- metanames
  for(i in 1:length(metanames)){
    if(metanames[i] %in% colnames(pd)) metaMatrix[,metanames[i]] <- pd[,metanames[i]]
  }
  
  if(all(is.na(metaMatrix[,"Sample name"]))) metaMatrix[,"Sample name"] <- sampleNames(normData)
 
  extraMeta <- pd[,setdiff(colnames(pd), metanames)]
  colnames(extraMeta) <- paste0("characteristics:",colnames(extraMeta))
  
  metaMatrix <- cbind(metaMatrix, extraMeta)
  
  newMeta <- metaTemplate
  newMeta[["SAMPLES"]] <- metaMatrix

  outfile <- paste0(basename, "meta.txt")
  
  cat("SERIES\n",file=outfile)
  sapply(newMeta["SERIES"], function(x) cat(paste0(as.character(x), "\n"), file=outfile,append=T))
  cat("SAMPLES\n",file=outfile,append=T)
  cat(paste0(paste(colnames(metaMatrix),collapse="\t"),"\n"),file=outfile,append=T)
  apply(metaMatrix, 1, function(x) cat(paste(x, collapse="\t"),"\n",file=outfile,append=T))
  cat("PROTOCOLS\n",file=outfile,append=T)
  sapply(newMeta["PROTOCOLS"], function(x) cat(paste0(as.character(x), "\n"), file=outfile,append=T))
  
  
  
}

createGEOMatrix <- function(normData,forceDetection=TRUE,basename="GEO",softTemplate=NULL,normalised=T){
  
  ###First create detection values if required.
  eMat <- exprs(normData)
  
  detMat <- matrix(nrow=nrow(eMat),ncol=ncol(eMat),NA)
  
  det <- Detection(normData)
  
  if(is.null(det)){
    
    if(forceDetection) {
      message("Calculating Detection scores")
      det <- try(calculateDetection(normData))
      
      if(class(det ) == "try-error") {
        message("Could not calculate detection scores")
      }
      else detMat <- det
    }
    
  } else detMat <- det
  
  Detection(normData) <- detMat
  
  
  
  
  if(is.null(softTemplate)){
  
  annoPkg <- paste("illumina", annotation(normData), ".db",sep="")
  
  
    annoLoaded <- require(annoPkg, character.only=TRUE)
    
      if(annoLoaded){
      
        mapEnv <-  as.name(paste("illumina", annotation(normData), "ARRAYADDRESS",sep=""))
      
        features <- mappedkeys(eval(mapEnv))
      
    
            
      }
  
    else message("Could not load required annotation package ",annoPkg)    
      
  }
  
  else{
    message("Parsing SOFT file....")
    
    if(length(grep(".gz", softTemplate) >0)) prbs <-unlist(lapply(strsplit(readLines(gzfile(softTemplate)), "\t"),function(x) x[1]))
    else prbs <-unlist(lapply(strsplit(readLines(softTemplate), "\t"),function(x) x[1]))
                   
    features <- prbs[grep("ILMN",prbs)] 
    message("Found ", length(features), " probes in SOFT file")
  }
  
    if(length(features) > length(featureNames(normData))){
      
      missingProbes <- setdiff(features, featureNames(normData))
      features <- intersect(features, featureNames(normData))
      message(length(missingProbes), " probes found in annotation, but not in normData: ", selectSome(missingProbes), ".See file MissingProbes.txt")
      write.table(missingProbes, file="MissingProbes.txt")
    } 
  
  extraProbes <- NA
  
  if(length(features) < length(featureNames(normData))){
    
    extraProbes <- setdiff(featureNames(normData),features)
    features <- intersect(features, featureNames(normData))
    message(length(extraProbes), " probes found in normData, but not in annotation: ", selectSome(extraProbes),".See file ExtraProbes.txt")
    write.table(fData(normData[extraProbes,]), file="ExtraProbes.txt")
  } 
  
    
  
      summaryData <- normData[features,]
    
      if(dim(summaryData)[1] == 0) message("Annotation map not successful")
      
      
   
      
      
      fullMat <- matrix(nrow=nrow(eMat), ncol=ncol(eMat)*2)
      fullMat[,seq(1, by = 2, len = ncol(eMat))] <- eMat
      fullMat[,seq(2, by = 2, len = ncol(eMat))] <- det
      colnames(fullMat)[seq(2, by = 2, len = ncol(eMat))] <- "Detection Pval"
      
      colnames(fullMat)[seq(1, by = 2, len = ncol(eMat))] <- colnames(eMat)
      fullMat <- cbind("ID_REF"= features, fullMat)
      
      outfile <- ifelse(normalised, paste0(basename, "ProcessedMatrix.txt"), paste0(basename, "RawMatrix.txt"))
      
      write.table(fullMat, file=outfile,sep="\t",quote=F,row.names=F)
      

  
}
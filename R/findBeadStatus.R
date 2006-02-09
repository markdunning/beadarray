findBeadStatus <- function(BLData, probes, array, log=FALSE, n=3, outputValid = FALSE, intProbeID = NULL){

  if(is.null(intProbeID)){
    intProbeID <- as.integer(BLData$ProbeID[,array])
  }
 
  outliers = valid = vector()
  
  for(i in 1:length(probes)){

    probe_ids = getProbeIndicesC(BLData, probe = probes[i], intProbe = intProbeID)
    inten <- BLData$R[probe_ids,array]
  
   #nas will be a list of beads which have NA intensity
    nas=NULL

    if(length(which(is.na(inten)))>0){

      nas = probe_ids[is.na(inten) | inten < 0]
      probe_ids = probe_ids[!is.na(inten)]
      inten = inten[!is.na(inten)]
    }

    if(log){
      raw_inten = log2(BLData$R[probe_ids,array])
      raw_inten = raw_inten[!is.na(raw_inten)]
    }
    else{
      raw_inten = inten
    }
  
    m = mean(inten, na.rm=TRUE, trim=0.5)
    ma = mad(inten, na.rm=TRUE)
  
  index = (inten > m + n *ma | inten < m - n*ma | raw_inten < 0)
  
    outliers.temp=probe_ids[index]

    outliers =c(outliers, outliers.temp, nas)
    if(outputValid)
      valid = c(valid, probe_ids[!index])
  }
  if(outputValid){
    result <- list(outliers = outliers, valid = valid)
    result
  }
  else
    outliers
}

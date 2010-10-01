
setGeneric("getArrayData", function(BLData, what="G", array=1, log=TRUE, method="illumina", n=3, trim=0.05)
   standardGeneric("getArrayData"))

setMethod("getArrayData", "BeadLevelList", function(BLData, what="G", array=1, log=TRUE, method="illumina", n=3, trim=0.05) {
   if(is.na(array))
      stop("'array' out of range")
   what = match.arg(what, choices=c("ProbeID", "GrnX", "GrnY", "G", "Gb", "R", "Rb", "wtsG", "wtsR", "residR", "residG", "M", "residM", "A", "beta", "wts", "xOffset", "yOffset"))

  if(what=="beta") {
     if(BLData@arrayInfo$channels=="two") {
       data=BLData[[array]][["R"]]/(BLData[[array]][["G"]]+BLData[[array]][["R"]])
     }
     else {
       stop("Need two-channel data to calculate per bead beta-values")
     }
   }
   if(what=="M") {
     if(BLData@arrayInfo$channels=="two") {
       data=log2.na(BLData[[array]][["R"]])-log2.na(BLData[[array]][["G"]])
     }
     else {
       stop("Need two-channel data to calculate per bead log-ratios")
     }
   }
   else if(what=="residM") {
     if(BLData@arrayInfo$channels=="two") 
       data = beadResids(BLData, what=gsub("resid", "", what), log=log, array=array, method=method, n=n, trim=trim)
      else
       stop("Need two-channel data to calculate per bead log-ratio residuals")
   }
   else if(what=="A") {
     if(BLData@arrayInfo$channels=="two")
       data=(log2.na(BLData[[array]][["R"]])+log2.na(BLData[[array]][["G"]]))/2
     else
       stop("Need two-channel data to calculate per bead average intensities")
   }
   else if(what=="residG") {
       data = beadResids(BLData, what=gsub("resid", "", what), log=log, array=array, method=method, n=n, trim=trim)
   }
   else if(what=="residR") {
     if(BLData@arrayInfo$channels=="two") 
       data = beadResids(BLData, what=gsub("resid", "", what), log=log, array=array, method=method, n=n, trim=trim)
      else
       stop("Need two-channel data to calculate per bead R residuals")
   }
   else if(what=="R" || what=="Rb") {
     if(BLData@arrayInfo$channels=="two") {
       data = BLData[[array]][[what]]
       if(log)
         data = log2.na(data)
     }
     else
       stop(paste("Need two-channel data to retrieve per bead", what, "values"))
   }
   else if(what =="wts"){
     data=BLData[[array]][["wts"]]
   }
   else if(what == "xOffset") {
     data=BLData@arrayInfo$xOffset[array]
   }
   else if(what == "yOffset") {
     data=BLData@arrayInfo$yOffset[array]
   }
else { # "G" or "Gb"
     data = BLData[[array]][[what]]
     if(log && (what=="G" || what=="Gb"))
       data = log2.na(data)
   }
   if(is.null(data))
     stop(paste("No", what, "data for array", array))  
   return(data)
 })
 
######################################################### 

setGeneric("arrayNames", function(object, arrays=NULL)
   standardGeneric("arrayNames"))

setMethod("arrayNames", "BeadLevelList", function(object, arrays=NULL) {
   if(is.null(arrays)) object@arrayInfo$arrayNames
   else object@arrayInfo$arrayNames[arrays]})

#########################################################

createBeadSummaryData = function(BLData, log = FALSE, imagesPerArray = 1, what="G", probes = NULL, arrays=NULL, method="illumina", n = 3, trim=0.05){
  arraynms = arrayNames(BLData)
  if((trim<0 || trim>0.5) && (method=="trim" || method=="winsorize"))
    stop("trim proportion must be between 0 and 0.5")
  if(method=="trim" && trim==0.5)
    method="median"
  if(!is.null(arrays) && !is.character(arrays))
    arraynms = arraynms[arrays]
  if(is.character(arrays))
    arraynms = which(arraynms %in% arrays)
  len = length(arraynms)
  what = match.arg(what, c("G", "R", "RG", "M", "A", "beta"))
  method = match.arg(method, c("illumina", "mean", "trim", "winsorize", "median"))
  method = match(method, c("illumina", "mean", "trim", "winsorize", "median"))
  whatelse = ""
  if(what=="RG") {
    if(BLData@arrayInfo$channels=="two") {
       what="G"
       whatelse="R"
     }
     else {
       stop("Need two-channel data to calculate summary R and G values")
     }
  }

  ##weights check (ensures backwards compatability)

  if(imagesPerArray == 1){
    sel = getArrayData(BLData, what="ProbeID", array=arraynms[1])!=0
    pr = getArrayData(BLData, what="ProbeID", array=arraynms[1])[sel]
    finten = getArrayData(BLData, what=what, log=log, array=arraynms[1])[sel]
    if("wts" %in% names(BLData[[arraynms[1]]]))
    {
      wts = getArrayData(BLData, what="wts", array=arraynms[1])[sel]
      pr = pr[wts==1]
      finten=finten[wts==1]     
    }
    nasinf = !is.finite(finten) | is.na(finten)
    pr = pr[!nasinf]
    finten = finten[!nasinf]
    binten = rep(0, length(finten))
  }
  else if(imagesPerArray == 2){
    if(length(arraynms)%%2!=0)
       stop("Need an even number of arrays when \'imagesPerArray=2\'")
    # order array by their pairs
    arrayord = order(arraynms)
    arraynms = arraynms[arrayord]
    # check that arrays are paired correctly
    tmp = unlist(strsplit(arraynms, "_"))
    chipnums = tmp[seq(1,length(tmp),by=3)]
    pos = tmp[seq(2,length(tmp),by=3)]
    stripnum = as.numeric(tmp[seq(3,length(tmp),by=3)])
    check = chipnums[seq(1,len, by=2)]==chipnums[seq(2,len, by=2)] & pos[seq(1,len, by=2)]==pos[seq(2,len, by=2)] & (stripnum[seq(1,len, by=2)]==stripnum[seq(2,len, by=2)]-1)
    if(sum(check)!=length(check))
       stop("Missing arrays")
    sel1 = getArrayData(BLData, what="ProbeID", array=arraynms[1])!=0
    sel2 = getArrayData(BLData, what="ProbeID", array=arraynms[2])!=0 
    pr = append(getArrayData(BLData, what="ProbeID", array=arraynms[1])[sel1],getArrayData(BLData, what="ProbeID", array=arraynms[2])[sel2])
    finten = append(getArrayData(BLData, what=what, log=log, array=arraynms[1])[sel1], getArrayData(BLData, what=what, log=log, array=arraynms[2])[sel2])

    if("wts" %in% names(BLData[[arraynms[1]]]) || "wts" %in% names(BLData[[arraynms[2]]]))
        {
      if("wts" %in% names(BLData[[arraynms[1]]])) {wts1 = getArrayData(BLData, what="wts", array=arraynms[1])[sel1]}
      else{wts1 = rep(1,length(sel1))}
      if("wts" %in% names(BLData[[arraynms[2]]])) {wts2 = getArrayData(BLData, what="wts", array=arraynms[2])[sel2]}
      else{wts2 = rep(1,length(sel2))}
      wts = append(wts1,wts2)
      #if(length(wts) != sum(wts) && method == 1){warning("Method = \"illumina\" used on an array with altered weights - have outliers been removed already? Consider using method = \"mean\" instead.")}
          pr = pr[wts==1]
          finten=finten[wts==1] 
        }
       nasinf = !is.finite(finten) | is.na(finten)

    pr = pr[!nasinf]
    finten = finten[!nasinf]
    binten = rep(0, length(finten))
    ord = order(pr)
    pr = pr[ord]
    finten = finten[ord]
  }
  else{
    stop("You can only specify 1 or 2 images per array")
  }
  
  if(is.null(probes)){
    probes = sort(unique(pr))
  }
    probes = probes[probes>0 & !is.na(probes)]
    noprobes = length(probes)

  if(imagesPerArray == 1){
    G  = GBeadStDev = GNoBeads = Gnooutliers = matrix(0,nrow = noprobes, ncol=len)
    colnames(G) = colnames(GBeadStDev) = colnames(GNoBeads) = colnames(Gnooutliers) = arraynms
    if(BLData@arrayInfo$channels=="two" && !is.null(BLData[[arraynms[1]]]$R) && whatelse=="R")
       R  = RBeadStDev = RNoBeads = Rnooutliers = G
    else R = NULL
  }

  else if(imagesPerArray == 2) {
    G =  GBeadStDev = GNoBeads = Gnooutliers = matrix(0,nrow = noprobes, ncol=(len/2))
    colnames(G) = colnames(GBeadStDev) = colnames(GNoBeads) = colnames(Gnooutliers) = arraynms[seq(1,len,by=2)]
    if(BLData@arrayInfo$channels=="two" && !is.null(BLData[[arraynms[1]]]$R) && whatelse=="R")
       R = RBeadStDev = RNoBeads = Rnooutliers = G
    else R = NULL
  }

  i = j = 1
   while(j <= len){

     probeIDs = as.integer(pr)

     start = 0
     blah = .C("createBeadSummary",  as.double(finten),  as.double(binten), probeIDs, as.integer(probes), as.integer(noprobes), as.integer(length(finten)),
                 fore = double(length = noprobes), back = double(length = noprobes), sd = double(length = noprobes), noBeads = integer(length = noprobes),
                 noOutliers = integer(length = noprobes), nextStart = as.integer(start), nmads = as.double(n), method = as.integer(method), trimfrac=as.double(trim), PACKAGE = "beadarray")

     G[,i] = blah$fore
     GBeadStDev[,i] = blah$sd
     GNoBeads[,i] = blah$noBeads
     Gnooutliers[,i] = blah$noOutliers
     
     if(BLData@arrayInfo$channels=="two" && !is.null(BLData[[arraynms[i]]]$R) && whatelse=="R") {
        if(imagesPerArray == 1){
           finten = getArrayData(BLData, what=whatelse, log=log, array=arraynms[i])[sel]
           nasinf = !is.finite(finten) | is.na(finten)
           finten = finten[!nasinf]
           binten = rep(0, length(finten))
        }
        else if(imagesPerArray == 2){
           finten = append(getArrayData(BLData, what=whatelse, log=log, array=arraynms[j])[sel1], getArrayData(BLData, what=whatelse, log=log, array=arraynms[j+1])[sel2])           
           nasinf = !is.finite(finten) | is.na(finten)
           finten = finten[!nasinf]
           binten = rep(0, length(finten)) 
           ord = order(pr)
           pr = pr[ord]
           finten = finten[ord]
        }

        start = 0
        blah = .C("createBeadSummary",  as.double(finten),  as.double(binten), probeIDs, as.integer(probes), as.integer(noprobes), as.integer(length(finten)),
                 fore = double(length = noprobes), back = double(length = noprobes), sd = double(length = noprobes), noBeads = integer(length = noprobes),
                 noOutliers = integer(length = noprobes), nextStart = as.integer(start), nmads = as.double(n), method = as.integer(method), trimfrac=as.double(trim), PACKAGE = "beadarray")

     R[,i] = blah$fore
     RBeadStDev[,i] = blah$sd
     RNoBeads[,i] = blah$noBeads
     Rnooutliers[,i] = blah$noOutliers
     }      
     j = j+imagesPerArray
     i = i + 1
     rm(probeIDs, blah)
     gc()
    
     if((imagesPerArray == 1) && (i <= len)) {    
       sel = getArrayData(BLData, what="ProbeID", array=arraynms[i])!=0
       pr = getArrayData(BLData, what="ProbeID", array=arraynms[i])[sel]   
       finten = getArrayData(BLData, what=what, log=log, array=arraynms[i])[sel]
       if("wts" %in% names(BLData[[arraynms[i]]]))
       {
         wts = getArrayData(BLData, what="wts", array=arraynms[i])[sel]
         pr = pr[wts==1]
         finten=finten[wts==1]  
       }
       nasinf = !is.finite(finten) | is.na(finten)
       pr = pr[!nasinf]
       finten = finten[!nasinf]
       binten = rep(0, length(finten))
     }
     else if((imagesPerArray == 2) && (j < len)) {
       sel1 = getArrayData(BLData, what="ProbeID", array=arraynms[j])!=0
       sel2 = getArrayData(BLData, what="ProbeID", array=arraynms[j+1])!=0 
       pr = append(getArrayData(BLData, what="ProbeID", array=arraynms[j])[sel1],getArrayData(BLData, what="ProbeID", array=arraynms[j+1])[sel2])
       finten = append(getArrayData(BLData, what=what, log=log, array=arraynms[j])[sel1], getArrayData(BLData, what=what, log=log, array=arraynms[j+1])[sel2])      
       if("wts" %in% names(BLData[[arraynms[j]]]) || "wts" %in% names(BLData[[arraynms[j+1]]]))
        {
          if("wts" %in% names(BLData[[arraynms[j]]])) {wts1 = getArrayData(BLData, what="wts", array=arraynms[j])[sel1]}
          else{wts1 = rep(1,length(sel1))}
          if("wts" %in% names(BLData[[arraynms[j+1]]])) {wts2 = getArrayData(BLData, what="wts", array=arraynms[j+1])[sel2]}
          else{wts2 = rep(1,length(sel2))}
          wts = append(wts1,wts2)

          pr = pr[wts==1]
              finten=finten[wts==1]     
            }
       nasinf = !is.finite(finten) | is.na(finten)
       pr = pr[!nasinf]
       finten = finten[!nasinf]
       binten = rep(0, length(finten))
       ord = order(pr)
       pr = pr[ord]
       finten = finten[ord]     
     }
   }
   
   
##Change the standard deviation to the standard error
if(method!=6) { # for all methods expect median
  GBeadStDev = GBeadStDev / sqrt(GNoBeads)
  if(!is.null(R))
    RBeadStDev = RBeadStDev / sqrt(RNoBeads)
}
  
if(whatelse=="R") { # create SNPSetIllumina
#    rownames(G) = rownames(R) = rownames(GBeadStDev) = rownames(RBeadStDev) = probes
#    require("beadarraySNP")
 #   BSData = new("SnpSetIllumina")
 #   assayData(BSData) = assayDataNew(G = G, R = R, storage.mode="list") # GBeadStDev = GBeadStDev, RBeadStDev = RBeadStDev,
        
        rownames(G) = rownames(R) = rownames(GBeadStDev) = rownames(RBeadStDev) = rownames(GNoBeads) = rownames(RNoBeads) = probes

        BSData = new("NChannelSet", R=R, G=G, GBeadStDev = GBeadStDev, RBeadStDev = RBeadStDev, GNoBeads=GNoBeads, RNoBeads=RNoBeads)  


}
else{
    BSData = new("ExpressionSetIllumina")
    assayData(BSData)=assayDataNew(exprs = G, se.exprs = GBeadStDev, NoBeads = GNoBeads,storage.mode="list")
    rownames(exprs(BSData)) = rownames(se.exprs(BSData)) = rownames(NoBeads(BSData)) = probes
    featureData(BSData) = new("AnnotatedDataFrame", data=data.frame(ProbeID=probes,row.names=probes))

}
  
if(nrow(pData(BLData))==len) {
  if(imagesPerArray==2)
    BSData@phenoData = new("AnnotatedDataFrame", data=pData(BLData@phenoData)[arrayord,,drop=FALSE][seq(1,len,by=2),,drop=FALSE])
  else
    BSData@phenoData = BLData@phenoData
  }
else {
  BSData@phenoData = new("AnnotatedDataFrame", data=data.frame(sampleName=colnames(G)))
}

  
if(!is.null(pData(BSData)$sampleName)) 
  sampleNames(BSData) = pData(BSData)$sampleName
else
  sampleNames(BSData) = colnames(G)


###Try and create pData properly for NChannelSet and make each column applicable to all channels


if(whatelse =="R"){
 varMetadata=data.frame(labelDescription=colnames(BSData@phenoData@data), channel="_ALL_")
 BSData@phenoData = new("AnnotatedDataFrame", data=data.frame(BSData@phenoData@data), varMetadata=varMetadata)
}

#if(length(BLData@annotation)==0) BSData@annotation="illuminaProbeIDs"
#else BSData@annotation=BLData@annotation
BSData@annotation=BLData@annotation

if("qcScores" %in% names(BLData@arrayInfo)) t=try(BSData@BeadLevelQC <- BLData@arrayInfo$qcScores,silent=TRUE)

BSData
} #)

# Define new bead level list class, which allows nBeads to differ between
# arrays, and stores the raw data in an environment to save on memory

setMethod("initialize", "BeadLevelList",
          function(.Object,  beadData = new.env(),
             phenoData=new("AnnotatedDataFrame", data=data.frame(name=I(character(0))),
               varMetadata=data.frame(labelDescription="Name in frame", row.names="name")),
             arrayInfo=list(arrayName=I(character(0)), nBeads=0), # probeindex=matrix(0),
             anno=character(0),
             beadAnno=data.frame(name=I(character(0))),
             scanMetrics=data.frame(name=I(character(0)))) {
               .Object@beadData = beadData
               .Object@phenoData = phenoData
               .Object@arrayInfo = arrayInfo
#               .Object@probeindex = probeindex
               .Object@annotation = anno
               .Object@beadAnno = beadAnno
               .Object@scanMetrics = scanMetrics
               .Object})


setMethod("[[",  # subsetting method, modified from the one in prada
  signature=c(x="BeadLevelList", i="ANY", j="missing"),
  definition=function(x, i, j="missing") {
    if(length(i)!=1)
      stop("subscript out of bounds (index must have length 1 in '[[')")
    if(is.numeric(i))
      i = arrayNames(x)[i]
    get(i, x@beadData, inherits=FALSE)
  },
  valueClass="data.frame")

setMethod("phenoData", "BeadLevelList", function(object) object@phenoData)

setMethod("pData", "BeadLevelList", function(object) pData(object@phenoData))

setGeneric("arrayNames", function(object, arrays=NULL)
   standardGeneric("arrayNames"))

setMethod("arrayNames", "BeadLevelList", function(object, arrays=NULL) {
   if(is.null(arrays)) object@arrayInfo$arrayNames
   else object@arrayInfo$arrayNames[arrays]})

setGeneric("numBeads", function(object, arrays=NULL)
   standardGeneric("numBeads"))

setMethod("numBeads", "BeadLevelList", function(object, arrays=NULL) {
   if(is.null(arrays)) object@arrayInfo$nBeads
   else object@arrayInfo$nBeads[arrays]}) 


setMethod("dim", "BeadLevelList", function(x) {
   arraynms = arrayNames(x)
   narrays = length(arraynms)
   nbeads = numBeads(x)
   c("nArrays"=narrays, "nBeads"=nbeads)
 })


setMethod("show", "BeadLevelList", function(object) {
   cat("Array information\n\n")
   show(object@arrayInfo)
   nrows = 5
   arraynms = arrayNames(object)
   cat(paste("Raw data from array", arraynms[1], "\n\n")) 
   show(object[[1]][1:nrows,])
   cat(paste("\n...", numBeads(object)[1]-nrows, "more rows of data\n\n"))
   cat(paste("... data for", length(arraynms)-1, "more array/s\n\n"))
   cat("Experimental information\n\n")
   show(pData(object@phenoData))
   cat("\nScanner metrics\n\n")
   show(summary(object@scanMetrics))
 })


setGeneric("getArrayData", function(BLData, what="G", array=1, log=TRUE, method="illumina", n=3, trim=0.05)
   standardGeneric("getArrayData"))

setMethod("getArrayData", "BeadLevelList", function(BLData, what="G", array=1, log=TRUE, method="illumina", n=3, trim=0.05) {
   if(is.na(array))
      stop("'array' out of range")
   what = match.arg(what, choices=c("ProbeID", "GrnX", "GrnY", "G", "Gb", "R", "Rb", "wtsG", "wtsR", "residR", "residG", "M", "residM", "A"))
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
   else { # "G" or "Gb"
     data = BLData[[array]][[what]]
     if(log && (what=="G" || what=="Gb"))
       data = log2.na(data)
   }
   if(is.null(data))
     stop(paste("No", what, "data for array", array))  
   return(data)
 })

setGeneric("copyBeadLevelList", function(object) 
  standardGeneric("copyBeadLevelList"))

setMethod("copyBeadLevelList", "BeadLevelList",
          function(object) {
             newobj = new("BeadLevelList")
             newobj@phenoData = object@phenoData
             newobj@arrayInfo = object@arrayInfo
#             newobj@probeindex = object@probeindex
             newobj@annotation = object@annotation
             newobj@beadAnno = object@beadAnno
             newobj@scanMetrics = object@scanMetrics
             arraynms = arrayNames(object)
             multiassign(arraynms, mget(arraynms, object@beadData),
                     envir = newobj@beadData, inherits=FALSE)
             newobj
          })


setGeneric("combineBeadLevelLists", function(object1, object2) 
  standardGeneric("combineBeadLevelLists"))

setMethod("combineBeadLevelLists", "BeadLevelList",
          function(object1, object2) {
             newbll = copyBeadLevelList(object1)
             newbll@arrayInfo$arrayNames = c(arrayNames(newbll), arrayNames(object2))
             newbll@arrayInfo$nBeads = c(numBeads(newbll), numBeads(object2))
             pData(newbll@phenoData) = rbind(pData(newbll@phenoData), pData(object2@phenoData))
#             newbll@probeindex = cbind(newbll1@probeindex , object2@probeindex)
             newbll@scanMetrics = rbind(newbll@scanMetrics, object2@scanMetrics)
             multiassign(arrayNames(object2), mget(arrayNames(object2), object2@beadData),
                     envir = newbll@beadData, inherits=FALSE)
             newbll
          })

#setGeneric("createBeadSummaryData", function(BLData, log = FALSE, n = 3, imagesPerArray = 2, probes = NULL)
#           standardGeneric("createBeadSummaryData"))

#setMethod("createBeadSummaryData", "BeadLevelList", function(BLData, log = FALSE, n = 3, imagesPerArray = 2, probes = NULL){

createBeadSummaryData = function(BLData, log = FALSE, imagesPerArray = 1, what="G", probes = NULL, arrays=NULL, method="illumina", n = 3, trim=0.05){
  arraynms = arrayNames(BLData)
  if((trim<0 || trim>0.5) && (method=="trim" || method=="winsorize"))
    stop("trim proportion must be between 0 and 0.5")
  if(!is.null(arrays) && !is.character(arrays))
    arraynms = arraynms[arrays]
  if(is.character(arrays))
    arraynms = which(arraynms %in% arrays)
  len = length(arraynms)
  what = match.arg(what, c("G", "R", "RG", "M", "A"))
  method = match.arg(method, c("illumina", "mean", "trim", "winsorize", "median"))
  method = match(method, c("illumina", "mean", "trim", "winsorize", "median"))
  if(method=="trim" && trim==0.5)
    method="median"
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
  if(imagesPerArray == 1){
    sel = getArrayData(BLData, what="ProbeID", array=arraynms[1])!=0
    pr = getArrayData(BLData, what="ProbeID", array=arraynms[1])[sel]
    finten = getArrayData(BLData, what=what, log=log, array=arraynms[1])[sel]
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
#    if(whatelse == "R") {
#       finten2 = append(getArrayData(BLData, what=whatelse, log=log, array=arraynms[1])[sel1], getArrayData(BLData, what=whatelse, log=log, array=arraynms[2])[sel2])
#       nasinf = !is.finite(finten) | is.na(finten) | !is.finite(finten2) | is.na(finten2)
#    }
#    else {
       nasinf = !is.finite(finten) | is.na(finten)
#    }
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
#     finten[!is.finite(finten) | is.na(finten)] = 0

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

#        finten[!is.finite(finten) | is.na(finten)] = 0

#        probeIDs = as.integer(pr)

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
       nasinf = !is.finite(finten) | is.na(finten)
       pr = pr[!nasinf]
       finten = finten[!nasinf]
       binten = rep(0, length(finten))
#       sel = BLData[[arraynms[i]]]$ProbeID!=0
#       pr = BLData[[arraynms[i]]]$ProbeID[sel]
#       finten = BLData[[arraynms[i]]]$G[sel]
#       binten = BLData[[arraynms[i]]]$Gb[sel]
     }
     else if((imagesPerArray == 2) && (j < len)) {
       sel1 = getArrayData(BLData, what="ProbeID", array=arraynms[j])!=0
       sel2 = getArrayData(BLData, what="ProbeID", array=arraynms[j+1])!=0 
       pr = append(getArrayData(BLData, what="ProbeID", array=arraynms[j])[sel1],getArrayData(BLData, what="ProbeID", array=arraynms[j+1])[sel2])
       finten = append(getArrayData(BLData, what=what, log=log, array=arraynms[j])[sel1], getArrayData(BLData, what=what, log=log, array=arraynms[j+1])[sel2])       
       nasinf = !is.finite(finten) | is.na(finten)
       pr = pr[!nasinf]
       finten = finten[!nasinf]
       binten = rep(0, length(finten))
       ord = order(pr)
       pr = pr[ord]
       finten = finten[ord]
#       binten = binten[ord]
#       sel1 = BLData[[arraynms[j]]]$ProbeID!=0
#       sel2 = BLData[[arraynms[j+1]]]$ProbeID!=0 
#       pr = append(BLData[[arraynms[j]]]$ProbeID[sel1], BLData[[arraynms[j+1]]]$ProbeID[sel2])
#       finten = append(BLData[[arraynms[j]]]$G[sel1], BLData[[arraynms[j+1]]]$G[sel2])
#       binten = append(BLData[[arraynms[j]]]$Gb[sel1], BLData[[arraynms[j+1]]]$Gb[sel2])
#       ord = order(pr)
#       pr = pr[ord]
#       finten = finten[ord]
#       binten = binten[ord]      
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
    require("beadarraySNP")
    BSData = new("SnpSetIllumina")
    assayData(BSData) = assayDataNew(G = G, R = R, storage.mode="list") # GBeadStDev = GBeadStDev, RBeadStDev = RBeadStDev,
    rownames(BSData@assayData[["G"]]) = rownames(BSData@assayData[["R"]]) = probes
#    rownames(BSData@assayData[["GBeadStDev"]]) = rownames(BSData@assayData[["RBeadStDev"]]) = probes

}
else{
    BSData = new("ExpressionSetIllumina")
    assayData(BSData)=assayDataNew(exprs = G, se.exprs = GBeadStDev, NoBeads = GNoBeads,storage.mode="list")
    rownames(exprs(BSData)) = rownames(se.exprs(BSData)) = rownames(NoBeads(BSData)) = probes

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

if(!is.null(BLData@annotation)) BSData@annotation="illuminaProbeIDs"
else BSData@annotation=BLData@annotation
BSData
} #)


#setGeneric("findAllOutliers", function(BLData, array=1, log = FALSE, n = 3, what="G") standardGeneric("findAllOutliers"))

#setMethod("findAllOutliers", "BeadLevelList", function(BLData, array=1, log = FALSE, n = 3, what="G"){
findAllOutliers = function(BLData, array=1, log = FALSE, n = 3, what="G"){
  probeList = getArrayData(BLData, array=array, what="ProbeID")
  probes = sort(unique(probeList[probeList>0]))
  inten = getArrayData(BLData, array=array, log=log, what=what)
  nasinf = is.na(inten) | !is.finite(inten)
  inten = inten[!nasinf]
  probeList = probeList[!nasinf]
  nbeads = length(inten)
  start = 0

  foo <- .C("findAllOutliers", as.double(inten), binStatus = integer(length = nbeads), as.integer(probeList), as.integer(probes), as.integer(length(probes)), as.integer(nbeads), as.integer(start), as.double(n), PACKAGE = "beadarray")

  which((probeList > 0) & (foo$binStatus == 0))
}# )

setGeneric("getProbeIntensities", function(BLData, ProbeIDs, array = 1, log = TRUE, what = "G") standardGeneric("getProbeIntensities"))

setMethod("getProbeIntensities", "BeadLevelList", function(BLData, ProbeIDs, array=1, log=TRUE, what="G"){
  return(getArrayData(BLData, array=array, what=what, log=log)[BLData[[array]]$ProbeID %in% ProbeIDs])
 })


findBeadStatus = function(BLData, probes, array=1, log=FALSE, what="G", n=3, outputValid = FALSE, intProbeID = NULL, ignoreList=NULL, probeIndex = NULL, startSearch = 1){

  if(is.null(intProbeID)){
    intProbeID = as.integer(sort(BLData[[array]]$ProbeID))
    probeIndex = c(1:length(intProbeID))
    probeIndex = probeIndex[sort.list(BLData[[array]]$ProbeID)]
  }

  outliers = valid = vector()

  for(i in 1:length(probes)){
    if(! (probes[i] %in% ignoreList) & sum(probes[i] %in% BLData[[array]]$ProbeID)>0){

      temp = getProbeIndicesC(BLData, probe = probes[i], array=array, intProbe = intProbeID, index = probeIndex, startSearch = startSearch)
      probe_ids = temp[[1]]
      startSearch = temp[[2]]

#      inten <- BLData@G[probe_ids,array]

   #nas will be a list of beads which have NA intensity
    nas=NULL

#    if(length(which(is.na(inten)))>0){

#      nas = probe_ids[is.na(inten) | inten < 0]
#   probe_ids = probe_ids[!is.na(inten)] 
#      inten = inten[!is.na(inten)]
#    }
      inten = getArrayData(BLData, array=array, what=what, log=log)[probe_ids]
#    if(log){
#      raw_inten = log2(BLData@G[probe_ids,array])
      inten = inten[!is.na(inten)]
      probe_ids = probe_ids[!is.na(inten)]
#    }
#    else{
#      raw_inten = inten
#    }

    m = mean(inten, na.rm=TRUE, trim=0.5) # this is the median
    ma = mad(inten, na.rm=TRUE) # apparently Illumina now use 3 SDs from the mean as default - check this
    index = (inten > m + n *ma | inten < m - n*ma | inten < 0)

    outliers.temp=probe_ids[index]

    outliers =c(outliers, outliers.temp, nas)
    if(outputValid)
      valid = c(valid, probe_ids[!index])
  }
  else cat("Probe with ID", probes[i], "ignored, or not found on array", array)
  }
  if(outputValid){
    result = list(outliers = outliers, valid = valid, nextStart = startSearch)
    result
  }
  else
    outliers
}


getProbeIndicesC = function(BLData, probe, array=1, intProbe, index, startSearch = 1){

  if(is.null(BLData[[array]]$ProbeID))
    stop("ProbeID column was not found in BeadLevelList object")
  nBeads = nrow(BLData[[array]])
  ind = .C("findIndices", as.integer(probe), intProbe, as.integer(nBeads), result = integer(length = 25000), pos = as.integer(startSearch), PACKAGE="beadarray")

  ind2 = vector()
  i = 1;
  while(ind$result[i] != 0){
    ind2 = c(ind2, index[ind$result[i]])
    i = i+1;
  }
  list(ind2, ind$pos)
}

# getUniqueProbeID <- function(BLData, array=1,...) {
#    sort(unique(BLData[[array]]$ProbeID),...)
#}


log2.na = function (x, ...)
{
    log2(ifelse(x > 0, x, NA), ...)
}

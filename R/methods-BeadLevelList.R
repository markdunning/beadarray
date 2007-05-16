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

setGeneric("arrayNames", function(object)
   standardGeneric("arrayNames"))

setMethod("arrayNames", "BeadLevelList", function(object) {
   object@arrayInfo$arrayNames})


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


setGeneric("getArrayData", function(BLData, which="G", array=1, log=TRUE)
   standardGeneric("getArrayData"))

setMethod("getArrayData", "BeadLevelList", function(BLData, which="G", array=1, log=TRUE) {
   which = match.arg(which, choices=c("ProbeID", "GrnX", "GrnY", "G", "Gb", "R", "Rb", "wtsG", "wtsR"))
   data = BLData[[array]][[which]]
   if(is.null(data))
     stop(paste("No", which, "data for array", array))
   if(log & (which=="G" | which=="Gb" | which=="R" | which=="Rb"))
     data = log2(data)
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
createBeadSummaryData = function(BLData, log = FALSE, n = 3, imagesPerArray = 2, probes = NULL){
arraynms = arrayNames(BLData)
  len = length(arraynms)

  if(imagesPerArray == 1){
    sel = BLData[[arraynms[1]]]$ProbeID!=0
    pr = BLData[[arraynms[1]]]$ProbeID[sel]
    finten = BLData[[arraynms[1]]]$G[sel]
    binten = BLData[[arraynms[1]]]$Gb[sel]
  }

  else if(imagesPerArray == 2){
    sel1 = BLData[[arraynms[1]]]$ProbeID!=0
    sel2 = BLData[[arraynms[2]]]$ProbeID!=0 
    pr = append(BLData[[arraynms[1]]]$ProbeID[sel1], BLData[[arraynms[2]]]$ProbeID[sel2])
    finten = append(BLData[[arraynms[1]]]$G[sel1], BLData[[arraynms[2]]]$G[sel2])
    binten = append(BLData[[arraynms[1]]]$Gb[sel1], BLData[[arraynms[2]]]$Gb[sel2])    
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
    if(BLData@arrayInfo$channels=="two" | !is.null(BLData[[arraynms[1]]]$R))
       R  = RBeadStDev = RNoBeads = Rnooutliers = G
    else R = NULL
  }

  else {
    G =  GBeadStDev = GNoBeads = Gnooutliers = matrix(0,nrow = noprobes, ncol=(len/2))
    colnames(G) = colnames(GBeadStDev) = colnames(GNoBeads) = colnames(Gnooutliers) = arraynms[seq(1,len,by=2)]
    if(BLData@arrayInfo$channels=="two" | !is.null(BLData[[arraynms[1]]]$R))
       R  = RBeadStDev = RNoBeads = Rnooutliers = G
    else R = NULL
  }

  i = j = 1
   while(j <= len){
    if(log){
     finten = log2(finten)
     binten = log2(binten)
     finten[!is.finite(finten) | is.na(finten)] = 0
     binten[!is.finite(binten) | is.na(binten)] = 0
   }

     probeIDs = as.integer(pr)

     start = 0
     blah = .C("createBeadSummary",  as.double(finten),  as.double(binten), probeIDs, as.integer(probes), as.integer(noprobes), as.integer(length(finten)),
                 fore = double(length = noprobes), back = double(length = noprobes), sd = double(length = noprobes), noBeads = integer(length = noprobes),
                 noOutliers = integer(length = noprobes), nextStart = as.integer(start), nmads = as.double(n), PACKAGE = "beadarray")

     G[,i] = blah$fore
     GBeadStDev[,i] = blah$sd
     GNoBeads[,i] = blah$noBeads
     Gnooutliers[,i] = blah$noOutliers
    
     if(BLData@arrayInfo$channels=="two" | !is.null(BLData[[arraynms[i]]]$R)) {
        if(imagesPerArray == 1){
           finten = BLData[[arraynms[i]]]$R[sel]
           binten = BLData[[arraynms[i]]]$Rb[sel]
        }
        else if(imagesPerArray == 2){
           finten = append(BLData[[arraynms[j]]]$R[sel1], BLData[[arraynms[j+1]]]$R[sel2])
           binten = append(BLData[[arraynms[j]]]$Rb[sel1], BLData[[arraynms[j+1]]]$Rb[sel2])    
  }
        if(log){
           finten = log2(finten)
           binten = log2(binten)
           finten[!is.finite(finten) | is.na(finten)] = 0
           binten[!is.finite(binten) | is.na(binten)] = 0
         }

        probeIDs = as.integer(pr)

        start = 0
        blah = .C("createBeadSummary",  as.double(finten),  as.double(binten), probeIDs, as.integer(probes), as.integer(noprobes), as.integer(length(finten)),
                 fore = double(length = noprobes), back = double(length = noprobes), sd = double(length = noprobes), noBeads = integer(length = noprobes),
                 noOutliers = integer(length = noprobes), nextStart = as.integer(start), nmads = as.double(n), PACKAGE = "beadarray")

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
       sel = BLData[[arraynms[i]]]$ProbeID!=0
       pr = BLData[[arraynms[i]]]$ProbeID[sel]
       finten = BLData[[arraynms[i]]]$G[sel]
       binten = BLData[[arraynms[i]]]$Gb[sel]
     }
     else if((imagesPerArray == 2) && (j < len)) {
       sel1 = BLData[[arraynms[j]]]$ProbeID!=0
       sel2 = BLData[[arraynms[j+1]]]$ProbeID!=0 
       pr = append(BLData[[arraynms[j]]]$ProbeID[sel1], BLData[[arraynms[j+1]]]$ProbeID[sel2])
       finten = append(BLData[[arraynms[j]]]$G[sel1], BLData[[arraynms[j+1]]]$G[sel2])
       binten = append(BLData[[arraynms[j]]]$Gb[sel1], BLData[[arraynms[j+1]]]$Gb[sel2])       
     }
   }

##Change the standard deviation to the standard error

GBeadStDev = GBeadStDev / sqrt(GNoBeads)
if(!is.null(R))
  RBeadStDev = RBeadStDev / sqrt(RNoBeads)

if(!is.null(R)) { # create SNPSetIllumina
#    rownames(G) = rownames(R) = rownames(GBeadStDev) = rownames(RBeadStDev) = probes
    BSData = new("SnpSetIllumina")
    assayData(BSData) = assayDataNew(G = G, R = R, GBeadStDev = GBeadStDev,
               RBeadStDev = RBeadStDev, storage.mode="list")
    rownames(BSData@assayData[["G"]]) = rownames(BSData@assayData[["R"]]) = probes = rownames(BSData@assayData[["GBeadStDev"]]) = rownames(BSData@assayData[["RBeadStDev"]]) = probes
    if(imagesPerArray==2)
      BSData@phenoData = new("AnnotatedDataFrame", data=pData(BLData@phenoData)[seq(1,len,by=2),,drop=FALSE])
    else
      BSData@phenoData = BLData@phenoData
}
else{
    BSData = new("ExpressionSetIllumina")
    assayData(BSData)=assayDataNew(exprs = G, BeadStDev=GBeadStDev, NoBeads = GNoBeads,storage.mode="list")
rownames(exprs(BSData)) = probes
    if(imagesPerArray==2)
      BSData@phenoData = new("AnnotatedDataFrame", data=pData(BLData@phenoData)[seq(1,len,by=2),,drop=FALSE])
    else
      BSData@phenoData = BLData@phenoData 
}
BSData@annotation=BLData@annotation
BSData
} #)


#setGeneric("findAllOutliers", function(BLData, array=1, log = FALSE, n = 3, which="G") standardGeneric("findAllOutliers"))

#setMethod("findAllOutliers", "BeadLevelList", function(BLData, array=1, log = FALSE, n = 3, which="G"){
findAllOutliers = function(BLData, array=1, log = FALSE, n = 3, which="G"){
  probes = sort(unique(BLData[[array]]$ProbeID[BLData[[array]]$ProbeID>0]))
  inten = getArrayData(BLData, array=array, log=log, which=which)
  inten[is.na(inten)| !is.finite(inten)] = 0

  probeList = BLData[[array]]$ProbeID
  nbeads = length(inten)
  start = 0

  foo <- .C("findAllOutliers", as.double(inten), binStatus = integer(length = nbeads), as.integer(probeList), as.integer(probes), as.integer(length(probes)), as.integer(nbeads), as.integer(start), as.double(n), PACKAGE = "beadarray")

  which((probeList > 0) & (foo$binStatus == 0))
}# )

setGeneric("getProbeIntensities", function(BLData, ProbeIDs, array = 1, log = TRUE, which = "G") standardGeneric("getProbeIntensities"))

setMethod("getProbeIntensities", "BeadLevelList", function(BLData, ProbeIDs, array=1, log=TRUE, which="G"){
  return(getArrayData(BLData, array=array, which=which, log=log)[BLData[[array]]$ProbeID %in% ProbeIDs])
 })


findBeadStatus = function(BLData, probes, array=1, log=FALSE, which="G", n=3, outputValid = FALSE, intProbeID = NULL, ignoreList=NULL, probeIndex = NULL, startSearch = 1){

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
      inten = getArrayData(BLData, array=array, which=which, log=log)[probe_ids]
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

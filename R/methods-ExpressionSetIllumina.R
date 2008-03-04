
setMethod("initialize", "ExpressionSetIllumina",
          function(.Object,
                   assayData = assayDataNew(exprs=exprs,se.exprs=se.exprs, NoBeads=NoBeads, Detection=Detection, Narrays=Narrays, arrayStDev=arrayStDev, DiffScore = DiffScore, storage.mode="list"),

                   phenoData = new("AnnotatedDataFrame"),
                   exprs=new("matrix"),
                   se.exprs=new("matrix"),
                   NoBeads=new("matrix"),
                   Detection=new("matrix"),
                   Narrays=new("matrix"),
                   arrayStDev =new("matrix"),
                   DiffScore=new("matrix"),
                   annotation = character(),
                   featureData = new("AnnotatedDataFrame"),
                   experimentData = new("MIAME"),

                   QCexprs = new("matrix"),
                   QCBeadStDev = new("matrix"),
                   QCNoBeads = new("matrix"),
                   controlType=new("matrix"),
                   QCData = assayDataNew(exprs = QCexprs,
                   se.exprs=QCBeadStDev, NoBeads=QCNoBeads,controlType=controlType,storage.mode="list"))
 {
            .Object<-callNextMethod(.Object,
                           assayData = assayData,
                           phenoData = phenoData,
                           experimentData = experimentData,
                           annotation = annotation,
                           featureData = featureData)
            .Object@QC=QCData
            .Object
          })


setMethod("[", "ExpressionSetIllumina", function(x, i, j, ..., drop = FALSE) {
          x<-callNextMethod(x, i, j, ..., drop=drop)
          if(!is.null(fData(x)) && !missing(i)) fData(x)<-fData(x)[i,, ..., drop = drop]
          x
})


setValidity("ExpressionSetIllumina", function(object) {
  assayDataValidMembers(assayData(object), c("exprs", "se.exprs", "NoBeads"))
})


setMethod("exprs", c("ExpressionSetIllumina"), function(object) assayDataElement(object, "exprs"))

setMethod("se.exprs", c("ExpressionSetIllumina"), function(object) assayDataElement(object, "se.exprs"))

setMethod("show", "ExpressionSetIllumina", function(object) {
  callNextMethod(object)
  
  cat("QC Information\n")
  cat(" Available Slots:  ")
  cat(names(object@QC))
  nms=selectSome(featureNames(object@QC))
  cat("\n  featureNames:", paste(nms, collapse=", "))
  nms=selectSome(sampleNames(object@QC))
  cat("\n  sampleNames:", paste(nms, collapse=", "))
  cat("\n")
})

setGeneric("QCInfo", function(object) standardGeneric("QCInfo"))

setMethod("QCInfo", "ExpressionSetIllumina", function(object) object@QC)

setGeneric("QC<-", function(object, value) standardGeneric("QC<-"))

setReplaceMethod("QC", "ExpressionSetIllumina", function(object, value){
	object@QC <- value
	object
})

setGeneric("NoBeads<-", function(object, value) standardGeneric("NoBeads<-"))

setReplaceMethod("NoBeads", "ExpressionSetIllumina", function(object, value){
	assayDataElementReplace(object, "NoBeads", value)
})

setGeneric("exprs<-", function(object, value) standardGeneric("exprs<-"))

setReplaceMethod("exprs", "ExpressionSetIllumina", function(object, value){
	assayDataElementReplace(object, "exprs", value)
})



setGeneric("NoBeads", function(object) standardGeneric("NoBeads"))

setMethod("NoBeads", "ExpressionSetIllumina", function(object) assayDataElement(object, "NoBeads"))

setGeneric("Detection<-", function(object, value) standardGeneric("Detection<-"))

setReplaceMethod("Detection", "ExpressionSetIllumina", function(object, value){
	assayDataElementReplace(object, "Detection", value)
})


setGeneric("Detection", function(object) standardGeneric("Detection"))

setMethod("Detection", "ExpressionSetIllumina", function(object) assayDataElement(object, "Detection"))

setGeneric("getVariance", function(object, offset=0) standardGeneric("getVariance"))

setMethod("getVariance", "ExpressionSetIllumina", function(object, offset=0) assayDataElement(object, "se.exprs")^2*assayDataElement(object, "NoBeads") + offset)


setReplaceMethod("exprs", c("ExpressionSetIllumina", "matrix"), function(object, value) {
  assayDataElementReplace(object, "exprs", value)
})

setReplaceMethod("se.exprs", c("ExpressionSetIllumina", "matrix"), function(object, value) {
  assayDataElementReplace(object, "se.exprs", value)
})


.mergeAssayData<-function(x, y, newdimnames) {
  # this is derived from assayData combine method
  # differences:
  # - allows different number of reporters/features
  # - will merge data from identical column names into 1 column ie rbind()) 
  # - only works on 2-dimensional assayData elements
  combineElement <- function(x, y) {
    outarr<-array(NA,dim=c(length(newdimnames[[1]]),length(newdimnames[[2]])),newdimnames)
    mode(outarr)<-mode(x)
    outarr[rownames(y),colnames(y)]<-y
    outarr[rownames(x),colnames(x)]<-x
    outarr
  }
  storage.mode <- storageMode(x)
  nmfunc <- assayDataElementNames

  if (storageMode(y) != storage.mode)
    stop(paste("assayData must have same storage, but are ",
               storage.mode, ", ", storageMode(y), sep=""))
  if (length(nmfunc(x)) != length(nmfunc(y)))
    stop("assayData have different numbers of elements:\n\t",
         paste(nmfunc(x), collapse=" "), "\n\t",
         paste(nmfunc(y), collapse=" "))
  if (!all(nmfunc(x) == nmfunc(y)))
    stop(paste("assayData have different element names:",
               paste(nmfunc(x), collapse=" "),
               paste(nmfunc(y), collapse=" "), sep="\n\t"))
               
  for (nm in nmfunc(x)) {
    x<-assayDataElementReplace(x,nm,combineElement(assayDataElement(x,nm),assayDataElement(y,nm)))
  }
  x
}

.mergePhenodata<-function(x , y, samples) {
  variables<-union(colnames(pData(x)),colnames(pData(y)))
  outarr<-array(data=NA,dim=c(length(samples),length(variables)),dimnames=list(samples,variables))
  outarr[sampleNames(y),colnames(pData(y))]<-as.matrix(pData(y))
  outarr[sampleNames(x),colnames(pData(x))]<-as.matrix(pData(x))
  pd<-data.frame(outarr)
  vardescs<-union(colnames(varMetadata(x)),colnames(varMetadata(y)))
  outarr<-array(data=NA,dim=c(length(variables),length(vardescs)),dimnames=list(variables,vardescs))
  outarr[colnames(pData(y)),colnames(varMetadata(y))]<-as.matrix(varMetadata(y))
  outarr[colnames(pData(x)),colnames(varMetadata(x))]<-as.matrix(varMetadata(x))
  vd<-data.frame(outarr)
  new("AnnotatedDataFrame", data=pd, varMetadata=vd)
}


setMethod("combine", c("ExpressionSetIllumina", "ExpressionSetIllumina"), function(x, y, ...) {
  if (class(x) != class(y))
    stop(paste("objects must be the same class, but are ",
               class(x), ", ", class(y), sep=""))
  ## we need a kind of merge functionality in order to combine OPA panels
  newdimnames<-list(union(featureNames(x),featureNames(y)),union(sampleNames(x),sampleNames(y)))
  x <- .mergeAssayData(x, y, newdimnames)
  # a bit of a hack to only keep the union, and discard double entries
  phenoData(x) <- .mergePhenodata(x, y, newdimnames[[2]])
  experimentData(x) <- combine(experimentData(x),experimentData(y))
    
  ## annotation -- constant
  if (any(annotation(x) != annotation(y))) {
    warning("objects have different annotations: ",
         annotation(x), ", ", annotation(y))
    annotation(x)<-unique(c(annotation(x),annotation(y)))
  }
  x
})

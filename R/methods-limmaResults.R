setMethod("initialize", "limmaResults",
          function(.Object,
                   assayData = assayDataNew(LogFC=LogFC,LogOdds=LogOdds, PValue=PValue, storage.mode="list"),
                   
                   phenoData = new("AnnotatedDataFrame"),
                   LogFC=new("matrix"),
                   LogOdds=new("matrix"),
                   PValue=new("matrix"),
                   annotation = character(),
                   featureData = new("AnnotatedDataFrame"),
                   experimentData = new("MIAME")
          )
{
            .Object<-callNextMethod(.Object,
                                    LogFC = LogFC,
                                    LogOdds = LogOdds,
                                    PValue = PValue,  
                                    experimentData = experimentData,
                                    annotation = annotation,
                                    featureData = featureData
            )
            
            .Object
          })


setGeneric("LogFC", function(object) standardGeneric("LogFC"))

setMethod("LogFC", signature(object="limmaResults"), function(object) assayDataElement(object, "LogFC"))

setGeneric("LogFC<-", function(object, value) standardGeneric("LogFC<-"))


setGeneric("LogOdds", function(object) standardGeneric("LogOdds"))

setMethod("LogOdds", signature(object="limmaResults"), function(object) assayDataElement(object, "LogOdds"))

setGeneric("LogOdds<-", function(object, value) standardGeneric("LogOdds<-"))



setGeneric("PValue", function(object) standardGeneric("PValue"))

setMethod("PValue", signature(object="limmaResults"), function(object) assayDataElement(object, "PValue"))

setGeneric("PValue<-", function(object, value) standardGeneric("PValue<-"))




setGeneric("ArrayWeights", function(object) standardGeneric("ArrayWeights"))

setMethod("ArrayWeights", signature(object="limmaResults"), function(object) object@ArrayWeights)

setGeneric("ArrayWeights<-", function(object, value) standardGeneric("ArrayWeights<-"))
setReplaceMethod("ArrayWeights",
                 signature=signature(
                   object="limmaResults",
                   value="numeric"),
                 function(object, value) {
                   object@ArrayWeights <- value
                   object
                 })




setGeneric("DesignMatrix", function(object) standardGeneric("DesignMatrix"))

setMethod("DesignMatrix", signature(object = "limmaResults"), function(object) object@DesignMatrix)
setGeneric("DesignMatrix<-", function(object, value) standardGeneric("DesignMatrix<-"))

setReplaceMethod("DesignMatrix",
                 signature=signature(
                   object="limmaResults",
                   value="matrix"),
                 function(object, value) {
                   object@DesignMatrix <- value
                   object
                 })


setGeneric("ContrastMatrix", function(object) standardGeneric("ContrastMatrix"))

setMethod("ContrastMatrix", signature(object = "limmaResults"), function(object) object@ContrastMatrix)
setGeneric("ContrastMatrix<-", function(object, value) standardGeneric("ContrastMatrix<-"))

setReplaceMethod("ContrastMatrix",
                 signature=signature(
                   object="limmaResults",
                   value="matrix"),
                 function(object, value) {
                   object@ContrastMatrix <- value
                   object
                 })

setMethod("show", signature(object="limmaResults"), function(object) {
  cat("Results of limma analysis\n")
  cat("Design Matrix used...\n")
  print(head(DesignMatrix(object)))
  cat(".....\n\n")
  cat("Array Weights....\n")
  cat("Contrast Matrix used...\n")
  print(head(ContrastMatrix(object)))
  cat("Array Weights....\n")
  cat(selectSome(round(ArrayWeights(object),digits=3)))
  cat("\nTop Table\n")
  
  topN <- order(LogOdds(object),decreasing=T)[1:10]
  sub <- object[topN,]
  df <- data.frame(fData(sub), LogFC=LogFC(sub), LogOdds=LogOdds(sub), pvalue = PValue(sub))
  print(df)
  
})






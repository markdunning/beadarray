limmaDE <- function(summaryData,SampleGroup,DesignMatrix=NULL,makeWts=TRUE,...){
  
  if(is.null(SampleGroup)) stop("You must define a SampleGroup for the differential expression\n")
  
  if (SampleGroup %in% colnames(pData(summaryData))) SampleGroup <- pData(summaryData)[,SampleGroup]
  else {
    print(paste(colnames(pData(summaryData)),collapse=" "))
    stop("The SampleGroup argument must be a column name in the phenoData slot. See above for list of valid strings")
  }
  
  design <- model.matrix(~0+as.factor(SampleGroup))
  colnames(design) <- as.character(levels(as.factor(SampleGroup)))
  
  
  
  
  contrast <- vector()
  
  for (a in 1:length(levels(as.factor(SampleGroup)))){
    for (b in 1:length(levels(as.factor(SampleGroup)))){ 
      if (a!=b){
        if (a<b){
          contrast[length(contrast)+1] <- paste(levels(as.factor(SampleGroup))[a],levels(as.factor(SampleGroup))[b],sep="-")
        }
      }
    }
  }
  
  if(makeWts){
    
    wts <- arrayWeights(exprs(summaryData),design);message("Calculating array weights")
    message("Array weights")
    wts
    
    fit <- lmFit(exprs(summaryData), design,weights=wts)
  }
  else   fit <- lmFit(exprs(summaryData), design)
  
  cnt <- paste(colnames(design)[1],colnames(design)[2],sep="-")
  cMat <- makeContrasts(contrasts =contrast,levels=design)
  
  
  fit2 <- contrasts.fit(fit,cMat)
  efit <- eBayes(fit2)
  
  
  deResults <- new("limmaResults",LogFC = efit$coef, PValue = efit$p.value, LogOdds = efit$lods,featureData=summaryData@featureData,phenoData = summaryData@phenoData)
  DesignMatrix(deResults) <- design
  ContrastMatrix(deResults) <- cMat
  ArrayWeights(deResults) <- wts
  annotation(deResults) <- annotation(summaryData)
  deResults@SampleGroup <- as.character(SampleGroup)
  deResults
  
}

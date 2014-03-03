

limmaDE <- function(summaryData,SampleGroup=NULL,DesignMatrix=NULL,makeWts=TRUE,...){
  
  if(is.null(SampleGroup)) SampleGroup <- SampleGroup(summaryData)
  
  
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
  
  
  deResults <- new("limmaResults",LogFC = efit$coef, PValue = efit$p.value, LogOdds = efit$lods,featureData=summaryData@featureData)
  DesignMatrix(deResults) <- design
  ContrastMatrix(deResults) <- cMat
  ArrayWeights(deResults) <- wts
  deResults
  
}







plotProbe <- function(symbol,rng, tx){
  data(genesymbol)  
  p1 <- autoplot(tx, which=genesymbol[symbol])
  p2 <- rng[which(rng %over% genesymbol[symbol])]
  tracks(p1, autoplot(p2,aes(fill=PROBEQUALITY)))
  
}




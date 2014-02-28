

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


setAs("ExpressionSet","ExpressionSetIllumina",
      function(from)
      {
        
        to <- new("ExpressionSetIllumina")
        exprs(to) <- exprs(from)
        phenoData(to) <- phenoData(from)  
        featureData(to) <- featureData(from)[,1:2]
        to@channelData[[1]] <- rep("G", length(sampleNames(to)))
        annotation(to) <- switch(annotation(from), 
                                 GPL6947="Humanv3", 
                                 GPL10558="Humanv4", 
                                 GPL6887="Mousev2", 
                                 GPL6102="Humanv2")
        to                          
      })


setAs("ExpressionSetIllumina", "GRanges",
    function(from)
    {
      
      
      annoName <- annotation(from)
      
      annoLoaded <- require(paste("illumina", annoName, ".db",sep=""), character.only=TRUE)
      
      if(annoLoaded){
        
        
        mapEnv <-  as.name(paste("illumina", annoName, "GENOMICLOCATION",sep=""))
        fn <- featureNames(from)
        fn <- fn[which(fn %in% mappedkeys(eval(mapEnv)))]
        
        locs <- mget(fn,eval(mapEnv),ifnotfound=NA)
        
        locs <- lapply(locs, function(x) gsub(" ", ",", x,fixed=T))
        
        asLocMatrix <- function(str){
          x<- do.call("rbind",sapply(strsplit(as.character(str), ",",fixed=T)[[1]], function(x) as.vector(strsplit(x, ":",fixed=T))))
        }
        
        locMat <- lapply(locs, asLocMatrix)
        
        rn <- rep(names(locs), unlist(lapply(locMat, nrow)))
        
        locMat <- do.call("rbind", locMat)
        rng <- GRanges(locMat[,1], IRanges(as.numeric(locMat[,2]), as.numeric(locMat[,3]),names=rn),strand=locMat[,4])
        #mcols(rng) <- df[match(names(rng), rownames(df)),]
        
        mcols(rng) <- data.frame(fData(from)[rn,], exprs(from)[rn,])

        sort(rng)
        
    }
      
    }

)



setAs("limmaResults", "GRanges",
      function(from)
      {
        annoName <- annotation(from)
        
        annoLoaded <- require(paste("illumina", annoName, ".db",sep=""), character.only=TRUE)
        
        if(annoLoaded){
          
          
          mapEnv <-  as.name(paste("illumina", annoName, "GENOMICLOCATION",sep=""))
          fn <- featureNames(from)
          fn <- fn[which(fn %in% mappedkeys(eval(mapEnv)))]
          
          locs <- mget(fn,eval(mapEnv),ifnotfound=NA)
          
          locs <- lapply(locs, function(x) gsub(" ", ",", x,fixed=T))
          
          asLocMatrix <- function(str){
            x<- do.call("rbind",sapply(strsplit(as.character(str), ",",fixed=T)[[1]], function(x) as.vector(strsplit(x, ":",fixed=T))))
          }
          
          locMat <- lapply(locs, asLocMatrix)
          
          rn <- rep(names(locs), unlist(lapply(locMat, nrow)))
          
          locMat <- do.call("rbind", locMat)
          rng <- GRanges(locMat[,1], IRanges(as.numeric(locMat[,2]), as.numeric(locMat[,3]),names=rn),strand=locMat[,4])
          #mcols(rng) <- df[match(names(rng), rownames(df)),]
          
          mcols(rng)$LogFC <- LogFC(from)[rn]
          mcols(rng)$LogOdds <- LogOdds(from)[rn]
          mcols(rng)$PValue <- PValue(from)[rn]
          sort(rng)
          
        }
        
      }
      
)


volcanoplot <- function(summaryData){
  
  data <- data.frame(LogFC = LogFC(summaryData), LogOdds = LogOdds(summaryData), fData(summaryData))
  colnames(data)[1] <- "LogFC"
  colnames(data)[2] <- "LogOdds"
  
  ggplot(data, aes(x = LogFC, y = LogOdds)) + geom_point(color="steelblue",alpha=0.3)
  
}

plotProbe <- function(symbol,rng, tx){
  data(genesymbol)  
  p1 <- autoplot(tx, which=genesymbol[symbol])
  p2 <- rng[which(rng %over% genesymbol[symbol])]
  tracks(p1, autoplot(p2,aes(fill=PROBEQUALITY)))
  
}




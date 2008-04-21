##Use normalize.invariant set to find a list of genes which do not change rank across arrays and then normalise to the 
##previously defined target distrbution


rankInvariantNormalise = function(exprs, T=NULL){

if(is.null(T)){

T = apply(exprs, 1, mean, na.rm=TRUE)

}


for (i in 1:ncol(exprs)){

curve = normalize.invariantset(exprs[,i], T)

exprs[,i] = predict(curve$n.curve, exprs[,i])$y
}

exprs

}





medianNormalise = function(exprs, log=TRUE){

exprs = as.matrix(exprs)

narrays = ncol(exprs)

if(log){
                          
exprs = log2(exprs)

}
                          
med = median(exprs,na.rm=TRUE)

for(i in 1:narrays){

exprs[,i]  = exprs[,i] - median(exprs[,i], na.rm=TRUE)

}


exprs = exprs + med


exprs

}


normaliseIllumina = function(BSData, method="quantile", transform="none", T=NULL, ...) {
  rownms = rownames(exprs(BSData))
  colnms = colnames(exprs(BSData))
  transform = match.arg(transform, c("none", "vst", "log2"))
  method = match.arg(method, c("quantile", "qspline", "vsn", "rankInvariant", "median", "none"))
  if(method=="vsn" && transform!="none"){
          cat(paste("\nmethod =", method, "not compatible with transform =", transform, "\nResetting transform = \"none\"\n\n"))
          transform="none"
  }
  if(transform=="log2") {
          BSData = assayDataElementReplace(BSData, "exprs", as.matrix(log2(exprs(BSData))))
  }
  else if(transform=="vst") {
          require("lumi")
          x = new("LumiBatch")
          x = assayDataElementReplace(x, "exprs", exprs(BSData))
          x = assayDataElementReplace(x, "se.exprs", se.exprs(BSData))
          x = assayDataElementReplace(x, "beadNum", NoBeads(BSData))
          if(!all(is.na(Detection(BSData))))
             x = assayDataElementReplace(x, "detection", Detection(BSData))
          BSData = assayDataElementReplace(BSData, "exprs", exprs(lumiT(x, method="vst", ...)))
          rm(x)
  }

  switch(method,
         quantile={
            require("affy")
            BSData = assayDataElementReplace(BSData, "exprs", normalize.quantiles(as.matrix(exprs(BSData))))
            rownames(exprs(BSData)) = rownms
            colnames(exprs(BSData)) = colnms
         },
         qspline={
            require("affy")
            BSData = assayDataElementReplace(BSData, "exprs", normalize.qspline(as.matrix(exprs(BSData))))
            rownames(exprs(BSData)) = rownms
            colnames(exprs(BSData)) = colnms
         },
         vsn={
            require("vsn")
            BSData = assayDataElementReplace(BSData, "exprs", exprs(vsn2(exprs(BSData))))
         },
         rankInvariant={
            BSData = assayDataElementReplace(BSData, "exprs", rankInvariantNormalise(exprs(BSData), T=T))

         },
         median={
            BSData = assayDataElementReplace(BSData, "exprs", medianNormalise(exprs(BSData), log=FALSE))
         })
  BSData
}

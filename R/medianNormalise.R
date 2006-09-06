"medianNormalise" <-
function(exprs, log=TRUE){

narrays = ncol(exprs)

exprs = log2(as.matrix(exprs))

med = median(exprs,na.rm=TRUE)

for(i in 1:narrays){

exprs[,i]  = exprs[,i] - median(exprs[,i], na.rm=TRUE)

}


exprs = exprs + med


exprs

}


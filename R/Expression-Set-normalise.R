##Use normalize.invariant set to find a list of genes which do not change rank across arrays and then normalise to the 
##previously defined target distrbution


rankInvariantNormalise = function(exprs, T=NULL){

if(is.null(T)){

T = apply(exprs, 1, mean)

}


for (i in 1:ncol(exprs)){

curve = normalize.invariantset(exprs[,i], T)

exprs[,i] = predict(curve$n.curve, exprs[,i])$y
}

exprs

}




"backgroundNormalise"<-function(object){

exprs = exprs(object)
QC = QCInfo(object)


if(length(which(colnames(QC$Signal) == "negative"))==0){

stop("Could not find the the column headings negative in the QC object")

}

col = which(colnames(QC$Signal) == "negative")

e = exprs - QC$Signal[,col]

e
 
}

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



##Use normalize.invariant set to find a list of genes which do not change rank across arrays and then normalise to the 
##previously defined target distrbution
rankInvariantNormalise = function(exprs, T){

for (i in 1:ncol(exprs)){

curve = normalize.invariantset(exprs[,i], T)

exprs[,i] = predict(curve$n.curve, exprs[,i])$y
}

exprs

}


backgroundNormalise = function(exprs, QC){

if(length(which(colnames(QC$Signal) == "negative"))==0){

stop("Could not find the the column headings negative in the QC object")

}

col = which(colnames(QC$Signal) == "negative")

e = exprs - QC$Signal[,col]

e

}

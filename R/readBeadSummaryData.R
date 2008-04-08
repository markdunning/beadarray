readBeadSummaryData = function(dataFile, qcFile=NULL, sampleSheet=NULL, header=TRUE, sep="\t", ProbeID="ProbeID",skip=8, columns = list(exprs = "AVG_Signal", se.exprs="BEAD_STDERR", NoBeads = "Avg_NBEADS", Detection="Detection Pval", Narrays="NARRAYS", arrayStDev = "ARRAY_STDEV"), qc.columns = list(exprs="AVG_Signal", se.exprs="BEAD_STDERR", NoBeads="Avg_NBEADS", Detection="Detection Pval", Narrays="NARRAYS", arrayStDev="ARRAY_STDEV"), controlID="ProbeID", annoPkg=NULL, 
qc.sep="\t", qc.skip=8, dec=".", quote="")
{

if(!(is.null(sampleSheet))){ 
samples = read.table(sampleSheet, sep=",", header=TRUE, skip=7)
}
  
r = read.table(as.character(dataFile), sep=sep, header=TRUE, skip=skip, dec=dec, quote=quote, as.is=TRUE, row.names=NULL, check.names=FALSE, strip.white=TRUE, comment.char="", fill=TRUE)

#foundColumns = NULL

index = grep(ProbeID, colnames(r))

if(length(index)!=0){

ProbeID = r[,index]
}


#check for non-unique probe names
if(length(ProbeID) != length(unique(ProbeID))){
      notdup = !duplicated(ProbeID)
      warning("ProbeIDs non-unique: consider setting 'ProbeID' to another column containing unique identifiers. ",   sum(!notdup), " repeated entries have been removed.\n")
      ProbeID = ProbeID[notdup]
      r = r[notdup,]
#    warning("ProbeIDs non-unique: consider setting 'ProbeID' to another column containing unique values.  Adding extension '.repX' to Xth replicate ID to enforce uniqueness.\n")
#    dups = unique(ProbeID[duplicated(ProbeID)])
#    for(j in 1:length(dups)) {
#      sel = ProbeID==dups[j]
#      ProbeID[sel] = paste(ProbeID[sel], ".rep", seq(1:sum(sel)), sep="")
#    }
}


data = list()

for(i in 1:length(columns)){


  index = grep(columns[[i]], colnames(r))

  if(length(index) == 0){
    cat("Could not find a column called: ", columns[[i]], "\n")
    data[[i]] = "" 
  }
  else{

#  foundColumns = append(foundColumns,names(columns)[i])  
  data[[i]] = r[,index]


  rownames(data[[i]]) = ProbeID
  colnames(data[[i]]) = sub(paste("(.|)", columns[[i]], "(.|)", sep=""), "", colnames(r)[index])
  if(!(is.null(sampleSheet))){
  colnames(data[[i]]) = as.character(samples[,4])
}
}
}

names(data) = names(columns) #foundColumns

BSData = new("ExpressionSetIllumina")

if(!is.null(annoPkg) && is.character(annoPkg))
  BSData@annotation = annoPkg

for(i in 1:length(data)){
  index = which(names(assayData(BSData))== names(data)[i])

  if(length(data[[i]])==1 && data[[i]]=="") {
    cat("No values stored in slot", names(data)[i], "\n")
  }

  else {
    assayData(BSData)[[index]] = as.matrix(data[[i]])
  }
}

if(!(is.null(qcFile))){
  BSData@QC = readQC(file=qcFile, sep=qc.sep, skip=qc.skip, columns=qc.columns, controlID=controlID, dec=dec, quote=quote)
  if(ncol(BSData@QC$exprs)!=ncol(exprs(BSData)))
     stop("Number of arrays doesn't agree: ", ncol(exprs(BSData)), " in main data set, versus ", ncol(BSData@QC$exprs), " in QC data.  Check your files.")
  notagree = colnames(BSData@QC$exprs)!=colnames(exprs(BSData))
  if(sum(notagree)!=0){
     for(i in 1:length(BSData@QC)) {
       reorder = sapply(colnames(BSData@QC[[i]]), FUN="grep", colnames(exprs(BSData)))
       BSData@QC[[i]] = BSData@QC[[i]][, reorder]
     }
  }
}

if(!(is.null(sampleSheet))){

p=new("AnnotatedDataFrame", samples,data.frame(labelDescription=colnames(samples), row.names=colnames(samples)))
phenoData(BSData) = p

}

BSData

}

readQC=function(file, columns=list(exprs="AVG_Signal", se.exprs="BEAD_STDERR", NoBeads="Avg_NBEADS", Detection="Detection Pval", Narrays="NARRAYS", arrayStDev="ARRAY_STDEV"), controlID = "ProbeID", sep="\t", skip=8, header=TRUE, dec=".", quote=""){

##
  
  r=read.table(as.character(file), sep=sep, header=header, skip=skip, quote=quote, as.is=TRUE, check.names=FALSE, strip.white=TRUE, comment.char="", fill=TRUE)


##If there is an ArrayID column, the QC file is BeadStudio version 1 and each row in the file is an array
#  read.table(file = fileName, header = TRUE, sep = sep, 
#        skip = nMetaDataLines, row.names = NULL, quote = "", 
#        as.is = TRUE, check.names = FALSE, strip.white = TRUE, 
#        comment.char = "", fill = TRUE)

  index = grep(controlID, colnames(r))
  
  if(length(index)!=0){
    ProbeID = r[,index]

    # check for non-unique probe names
    if(length(ProbeID) != length(unique(ProbeID))){
      notdup = !duplicated(ProbeID)
      warning("controlIDs non-unique: ", sum(!notdup), " repeated entries have been removed.\n")
      ProbeID = ProbeID[notdup]
      r = r[notdup,]
    }
  }
  else { # can't find controlIDs - report warning
    warning("controlID does not have a match in", file, "check, qc.skip and controlID arguments.  Using numbers as rownames")
    ProbeID = seq(1:nrow(r))
  }
     
 ArrayID = grep("ArrayID", colnames(r))
 
 if(length(ArrayID) != 0){
   
 BeadStudioVersion=1 
  
}


 else BeadStudioVersion=2 

  
#  foundColumns = NULL


data = list()

# count = 1
for(i in 1:length(columns)){

  index = grep(columns[[i]], colnames(r))

  if(length(index) == 0){
    cat("[readQC] Could not find a column called: ", columns[[i]], "\n")
    data[[i]] = ""
  }
  else{
#  foundColumns = append(foundColumns,names(columns)[i], "\n")  

  if(BeadStudioVersion == 2) { 
    data[[i]] = as.matrix(r[,index])
  }
  else {
    data[[i]] = as.matrix(t(r[,index]))
    colnames(data[[i]]) = r[,ArrayID]
  }           
#  count = count+1
  rownames(data[[i]]) = ProbeID
  colnames(data[[i]]) = sub(paste("(.|)", columns[[i]], "(.|)", sep=""), "", colnames(r)[index])
}
}

names(data) = names(columns) # foundColumns
QC = assayDataNew(exprs=new("matrix"), se.exprs=new("matrix"), Detection=new("matrix"), NoBeads=new("matrix"), storage.mode="list")
  

for(i in 1:length(data)){

  if(length(data[[i]])==1 && data[[i]]=="") {
    cat("[readQC] No values stored in slot", names(data)[i], "\n")
  }
  else {
    index = which(names(QC)== names(data)[i])
    if(length(index)>0){
      QC[[index]] = data[[i]]
    }
  }

}

 QC

}

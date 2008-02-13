readBeadSummaryData<- function(dataFile, qcFile=NULL, sampleSheet=NULL, header=TRUE, sep="\t", ProbeID="ProbeID",skip=8, columns = list(exprs = "AVG_Signal", NoBeads = "Avg_NBEADS", Detection="Detection", se.exprs="BEAD_STDERR", Narrays="NARRAYS", arrayStDev = "ARRAY_STDEV"), qc.columns = list(controlID="ProbeID", controlType="TargetID", exprs="AVG_Signal", Detection="Detection", Narrays="NARRAYS", se.exprs="BEAD_STDERR", NoBeads="Avg_NBEADS", arrayStDev="ARRAY_STDEV"), annoPkg=NULL, qc.sep="\t", qc.skip=8,...)
{

if(!(is.null(sampleSheet))){ 
samples = read.table(sampleSheet, sep=",", header=T, skip=7)
}
  
r=read.table(as.character(dataFile), sep=sep, header=T, skip=skip,...) # quote="")


foundColumns = NULL

index = grep(ProbeID, colnames(r))

if(length(index)!=0){

ProbeID = r[,index]
}


data = list()


for(i in 1:length(columns)){


  index = grep(columns[[i]], colnames(r))

  if(length(index) == 0){

    cat(paste("Could not find a column called: ", columns[[i]]))
    cat("\n")
  }
  else{

  foundColumns = append(foundColumns,names(columns)[i])  
  data[[i]] = r[,index]

  ##check for unique probe names
  if(length(ProbeID) == length(unique(ProbeID))){

  rownames(data[[i]]) = ProbeID
}

  if(!(is.null(sampleSheet))){
  
  colnames(data[[i]]) = as.character(samples[,4])
}
}
}

names(data) = foundColumns



BSData = new("ExpressionSetIllumina")
if(!is.null(annoPkg) && is.character(annoPkg))
  BSData@annotation = annoPkg
for(i in 1:length(data)){

  index = which(names(assayData(BSData))== names(data)[i])

  if(length(index)>0){
  assayData(BSData)[[index]] = data[[i]]
}
    else{
    cat(paste("Did not find a slot called :", names(data)[i]))
    cat("\n")
  }

}


if(!(is.null(qcFile))){

BSData@QC = readQC(file=qcFile, sep=qc.sep, skip=qc.skip, columns=qc.columns)


}


if(!(is.null(sampleSheet))){

p=new("AnnotatedDataFrame", samples,data.frame(labelDescription=colnames(samples), row.names=colnames(samples)))
phenoData(BSData) = p

}

BSData


}

readQC=function(file, columns = list(exprs = "AVG_Signal", NoBeads = "Avg_NBEADS", Detection="Detection", se.exprs="BEAD_STDERR", Narrays="NARRAYS", arrayStDev = "ARRAY_STDEV", controlID = "ProbeID", controlType="TargetID"),sep="\t",skip=7,header=T){


##
  
  r=read.table(as.character(file), sep=sep, header=header, skip=skip)


##If there is an ArrayID column, the QC file is BeadStudio version 1 and each row in the file is an array
#  read.table(file = fileName, header = TRUE, sep = sep, 
#        skip = nMetaDataLines, row.names = NULL, quote = "", 
#        as.is = TRUE, check.names = FALSE, strip.white = TRUE, 
#        comment.char = "", fill = TRUE)
 

  
 ArrayID = grep("ArrayID", colnames(r))

 if(length(ArrayID) != 0){
   
 BeadStudioVersion=1 
  
}


 else BeadStudioVersion=2 

  
  foundColumns = NULL


data = list()


for(i in 1:length(columns)){


  index = grep(columns[[i]], colnames(r))

  if(length(index) == 0){

    cat(paste("Could not find a column called: ", columns[[i]]))
    cat("\n")
  }
  else{

  foundColumns = append(foundColumns,names(columns)[i])  

  if(BeadStudioVersion == 2) data[[i]] = as.data.frame(r[,index])
  else{
    data[[i]] = as.data.frame(t(r[,index]))
    colnames(data[[i]]) = r[,ArrayID]
  }           
  


}
}

names(data) = foundColumns
QC = assayDataNew(exprs=new("matrix"), se.exprs=new("matrix"), Detection=new("matrix"), NoBeads=new("matrix"),controlType=new("matrix"),storage.mode="list")
  

for(i in 1:length(data)){

  index = which(names(QC)== names(data)[i])

  if(length(index)>0){
  QC[[index]] = data[[i]]
}
    else{
    cat(paste("Did not find a slot called :", names(data)[i]))
    cat("\n")
  }
}

 QC
  

}

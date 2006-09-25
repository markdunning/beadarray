
readQC=function(file, columns=list(Biotin="AVG.Signal.biotin", cy3_high="AVG.Signal.cy3_hyb_high", cy3_low="AVG.Signal.cy3_hyb_low", cy3_med="AVG.Signal.cy3_hyb_med", gene="AVG.Signal.gene", hs="AVG.Signal.high_stringency_hyb", house="AVG.Signal.housekeeping", labeling="AVG.Signal.labeling", mm="AVG.Signal.low_stringency_hyb_mm2", pm="AVG.Signal.low_stringency_hyb_pm", negative="AVG.Signal.negative"),skip=7,sep=",",header=T){


  r=read.table(as.character(file), sep=sep, header=header, skip=skip)


   signal = matrix(nrow=nrow(r), ncol=length(columns))

  for(i in 1:length(columns)){

    m = which(as.character(colnames(r))==columns[i])

    signal[,i] = r[,m]

  }

  rownames(signal)  = as.character(r[,1])

  colnames(signal) = as.character(names(columns))



  var = matrix(nrow=nrow(r), ncol=length(columns))

  for(i in 1:length(columns)){

    m = which(as.character(colnames(r))==columns[i]) + 1

    var[,i] = r[,m]

  }

  rownames(var)  = as.character(r[,1])

  colnames(var) = as.character(names(columns))



  detection = matrix(nrow=nrow(r), ncol=length(columns))

  for(i in 1:length(columns)){

    m = which(as.character(colnames(r))==columns[i]) + 2

    detection[,i] = r[,m]

  }

  rownames(detection)  = as.character(r[,1])

  colnames(detection) = as.character(names(columns))

  QC = assayDataNew(Signal = signal, StDev=var, Detection=detection, storageMode="list")

  QC
}

"plotQC"<-function(object){


  if(class(BSData)== "ExpressionSetIllumina") QC = QCInfo(BSData)
  else if(class(BSData) != "AssayData") stop("Input object must of type AssayData or ExpressionSetIllumina")

  
  par(mfrow=c(3,2))

  plotQC.figure1(QC,...)
  plotQC.figure2(QC,...)
  plotQC.figure3(QC,...)
  plotQC.figure4(QC,...)
  plotQC.figure5(QC,...)
  plotQC.figure6(QC,...)

}




plotQC.figure1 = function(QC,log=FALSE,...){
  if(class(BSData)== "ExpressionSetIllumina") QC = QCInfo(BSData)
  else if(class(BSData) != "AssayData") stop("Input object must of type AssayData or ExpressionSetIllumina")

  low = which(colnames(QC$Signal) == "cy3_low")
  med = which(colnames(QC$Signal) == "cy3_med")
  high = which(colnames(QC$Signal) == "cy3_high")

  boxplot(as.data.frame(QC$Signal[,c(low,med,high)]), names=c("low","med","high"), ylab="Signal", main="Hybridisation Controls",...)

 }

plotQC.figure2 = function(QC,log=FALSE,...){

  if(class(BSData)== "ExpressionSetIllumina") QC = QCInfo(BSData)
  else if(class(BSData) != "AssayData") stop("Input object must of type AssayData or ExpressionSetIllumina")
   
  
  bac = which(colnames(QC$Signal) == "negative")

 boxplot(QC$Signal[,bac], main="Background",...)
  
}


plotQC.figure3 = function(QC, log=FALSE,...){

  if(class(BSData)== "ExpressionSetIllumina") QC = QCInfo(BSData)
  else if(class(BSData) != "AssayData") stop("Input object must of type AssayData or ExpressionSetIllumina")
   

  biotin = which(colnames(QC$Signal) == "Biotin")
  hs = which(colnames(QC$Signal) == "hs")

  boxplot(as.data.frame(QC$Signal[,c(biotin, hs)]), names=c("biotin", "high stringency"), ylab="Signal", main="Biotin vs High Stringency",...)

}


plotQC.figure4 = function(QC, log=FALSE,...){

  if(class(BSData)== "ExpressionSetIllumina") QC = QCInfo(BSData)
  else if(class(BSData) != "AssayData") stop("Input object must of type AssayData or ExpressionSetIllumina")
   

  hkp = which(colnames(QC$Signal) == "house")
  gene = which(colnames(QC$Signal) == "gene")

  boxplot(as.data.frame(QC$Signal[,c(hkp, gene)]), names=c("housekeeping", "gene"), ylab="Signal", main="Housekeeping vs All Genes",...)
}

plotQC.figure5 = function(QC, log=FALSE,...){

  if(class(BSData)== "ExpressionSetIllumina") QC = QCInfo(BSData)
  else if(class(BSData) != "AssayData") stop("Input object must of type AssayData or ExpressionSetIllumina")
   

  mm = which(colnames(QC$Signal) == "mm")
  pm = which(colnames(QC$Signal) == "pm")

  boxplot(as.data.frame(QC$Signal[,c(mm,pm)]), names=c("mm2", "pm"), ylab="Signal", main="Low Stringency",...)

}

plotQC.figure6 = function(QC, log=FALSE,...){

  if(class(BSData)== "ExpressionSetIllumina") QC = QCInfo(BSData)
  else if(class(BSData) != "AssayData") stop("Input object must of type AssayData or ExpressionSetIllumina")
   

  lab = which(colnames(QC$Signal) == "labeling")
  bac = which(colnames(QC$Signal) == "negative")

  boxplot(as.data.frame(QC$Signal[,c(lab, bac)]), names=c("labeling", "background"), main="Labeling vs Background",...)
}


           
           
"singleQCPlot"<-function(object, type="Biotin", log=FALSE, whatToPlot="Signal"){

  if(class(object)== "ExpressionSetIllumina") QC = QCInfo(BSData)
  else if(class(BSData) != "AssayData") stop("Input object must of type AssayData or ExpressionSetIllumina")
   

  if(whatToPlot == "Signal"){
  data = QC$Signal
  col = which(colnames(QC$Signal)== type)
  }
else if(whatToPlot == "Var"){
  data=QC$Var
  col = which(colnames(QC$Var)== type)   	
}	
else if(whatToPlot == "Detection"){
  data=QC$Detection
  col = which(colnames(QC$Detection)== type)   	
}



  if(length(col) == 0){

    stop(paste("type argument must be one of: ", colnames(QC$Signal)))
  }
  if(log){
  plot(log2(data[,col]), xlab="ArrayIndex", ylab="Signal",...)
  lines(log2(data[,col]))
  }

else{
   plot(data[,col], xlab="ArrayIndex", ylab="Signal",...)
   lines(data[,col])
 }

}


  


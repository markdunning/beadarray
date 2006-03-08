
plotProbeVariation = function(BSData, id, log=TRUE, ylim=range(0,16),...){
 

  values = as.matrix(BSData$R[BSData$ProbeID==id,])

  if(log) values=log2(values)
  
  plot(values)


}

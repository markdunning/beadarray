
plotProbeVariation = function(BSData, ProbeID, log=TRUE, ylim=range(0,16), main=ProbeID,...){


  values = BSData$R[BSData$ProbeID[,1]==ProbeID,]

  if(log) values=log2(values)
  
  plot(values, type="l", ylim=ylim, main=main,...)


}

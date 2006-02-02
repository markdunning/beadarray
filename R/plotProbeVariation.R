
plotProbeVariation = function(BSData, probeID, log=TRUE, ylim=range(0,16), main=probeID,...){


  values = BSData$R[BSData$probeID[,1]==probeID,]

  if(log) values=log2(values)
  
  plot(values, type="l", ylim=ylim, main=main,...)


}

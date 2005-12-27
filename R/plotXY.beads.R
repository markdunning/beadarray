"plotXY.beads" <-
function(BSData,array1, array2=0, log=FALSE, identify=FALSE,label=FALSE, fold=2, highlight=NULL){

#XY plot of either two samples against each other, or red and green channels of one channel


  if(!class(BSData) == "BeadSummaryList"){
    stop("BeadSummaryList object required!")
  }  

  probes = sort(unique(BSData$probeID[BSData$probeID[,1] > 0,1]))

  if (array2!=0){

    if(log){

      x = log2(BSData$R[,array1])
      y = log2(BSData$R[,array2])


      xmax = 20
      xbox=18
      yspacing=0.3
    }

    else{
      x = BSData$R[,array1]
      y = BSData$R[,array2]

      xmax=2^17
      xbox=100000
      yspacing=3000
    }
  }
  else{
    if(log){
      x = log2(BSData$R[,array1])
      y = log2(BSData$G[,array1])

      xmax = 20
      xbox=18
      yspacing=0.3
    }
    else{
      x = BSData$R[,array1]
      y = BSData$G[,array1]

      xmax=2^17
      xbox=100000
      yspacing=3000
    }
  }

  plot(x,y, col="black",xlim=range(min(x),xmax))

  abline(0,1)

  if(label){
    
    for(i in 1:nrow(spottypes)){

      x1 = x[BSData$probeID[,1]==spottypes[i,2]]
      y1 = y[BSData$probeID[,1]==spottypes[i,2]]

      points(x1,y1, col=as.character(spottypes[i,4]))

      types = as.character(unique(spottypes[,1]))

      cols = as.character(unique(spottypes[,4]))

      legend(min(x),max(y), types,cols,cex=0.8)
    }
  }


if(log){



abline(log2(fold),1)
abline(-log2(fold),1)

}

else{
abline(0,fold)
abline(0,1/fold)


}

if(length(highlight!=0)){

for(i in 1:length(highlight)){

x1 = x[BSData$probeID[,1]==highlight[i]]
y1 = y[BSData$probeID[,1]==highlight[i]]


points(x1,y1, cex=1.6, col="red", pch=3)

}
}

if(identify){

id=identify(x,y, x,n=1, plot=FALSE)

info=getProbeInfo(BSData$genes, probes[id])

if(nrow(info) != 0){ 

print(info)
text(xbox,max(y)-yspacing,paste("Probe ID: ", probes[id],sep=""))
text(xbox,max(y)-(yspacing*2),paste("Gene ID: ", as.character(info[,2]), sep=""))
text(xbox,max(y)-(yspacing*3),paste("Search Key:", as.character(info[,1]),sep=""))
text(xbox,max(y)-(yspacing*4),paste("Value from SAM 1: ", round(x[id],2), sep=""))
text(xbox,max(y)-(yspacing*5),paste("Value from SAM 2: ", round(y[id],2), sep=""))

}

else{

info=spottypes[spottypes$ID==probes[id],]

print(info)

text(xbox,max(y)-yspacing,paste("Probe ID: ", probes[id],sep=""))
text(xbox,max(y)-(yspacing*2),paste("Control Type: ", as.character(info[,3]),sep=""))
text(xbox,max(y)-(yspacing*3),paste("Source Reference: ", as.character(info[,2]),sep=""))
text(xbox,max(y)-(yspacing*4),paste("Value from SAM 1: ", round(x[id],2), sep=""))
text(xbox,max(y)-(yspacing*5),paste("Value from SAM 2: ", round(y[id],2), sep=""))

}



}



}


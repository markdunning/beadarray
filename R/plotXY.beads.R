"plotXY.beads" <-
function(exprs,array1, array2=0, log=FALSE, identify=FALSE,label=FALSE, fold=2, highlight=NULL, sampleSize=NULL,...){

#XY plot of either two samples against each other, or red and green channels of one channel

exprs=as.matrix(exprs)

  if (array2!=0){

    if(log){

      x = log2(exprs[,array1])
      y = log2(exprs[,array2])


      xmax = 16
      xbox=18
      yspacing=0.3
    }

    else{
      x = exprs[,array1]
      y = exprs[,array2]

      xmax=2^17
      xbox=100000
      yspacing=3000
    }
  }
  else{
    if(log){
      x = log2(exprs[,array1])
      y = log2(exprs[,array1])

      xmax = 16
      xbox=18
      yspacing=0.3
    }
    else{
      x = exprs[,array1]
      y = exprs[,array1]

      xmax=2^17
      xbox=100000
      yspacing=3000
    }
  }
 if(!is.null(sampleSize)){
 s = sample(1:length(x), sampleSize)
 x=x[s]
 y=y[s]
 }
	

  plot(x,y, col="black",xlim=range((max(0,min(x),na.rm=TRUE)),16), xlab = "", ylab = "", pch = 16, cex = 0.4, ...)

  abline(0,1)

  if(label){

    status = BSData$genes$Status


    values <- attr(status, "values")

    nvalues = length(values)

    
    sel <- !(status %in% values)

    col <- attr(status, "Colour")

    col <- rep(col, length=nvalues)

    for(i in 1:nvalues){
      sel <- status==values[i]    

       ids = BSData$genes$ProbeID[sel]
      
      points(x[which(BSData$ProbeID %in% ids)], y[which(BSData$ProbeID %in% ids)], col=col[i])

    }

    
    
  }


if(log){

abline(log2(fold),1, lty=2)
abline(-log2(fold),1, lty=2)

}

else{
abline(log2(fold),1, lty=2)
abline(-log2(fold),1, lty=2)


}

if(length(highlight!=0)){

for(i in 1:length(highlight)){

x1 = x[BSData$ProbeID==highlight[i]]
y1 = y[BSData$ProbeID==highlight[i]]


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


"plotMA.beads" <-
function(BSData, array1, array2=0, identify=FALSE, label=FALSE, highlight=NULL, log=TRUE, main=NULL, ma.ylim=2,...){

  if(!class(BSData) == "BeadSummaryList"){
    stop("BeadSummaryList object required!")
  }

probes = sort(unique(BSData$ProbeID[BSData$ProbeID > 0]))

if(array2!=0){

if(log){

x = 0.5*(log2(BSData$R[,array1]) + log2(BSData$R[,array2]))
y = log2(BSData$R[,array1])- log2(BSData$R[,array2])
  

}

else{

x = 0.5*(BSData$R[,array1] + BSData$R[,array2])
y = BSData$R[,array1]- BSData$R[,array2]
}


}

else{
x = log2(BSData$R[,array1])
y = log2(BSData$G[,array1])

}

  plot(x,y, pch=16,cex=0.4, ylim=range(ma.ylim,-ma.ylim),xlim=range(5,16), main=main, xlab = "", ylab = "", ...) 

  abline(h=c(-1,0,1),lty=c(2,1,2)) 
  
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

if(length(highlight!=0)){

for(i in 1:length(highlight)){

x1 = x[BSData$ProbeID[,1]==highlight[i]]
y1 = y[BSData$ProbeID[,1]==highlight[i]]


points(0.5*(x1+y1),y1-x1, cex=1.6, col="red", pch=3)

}


}



if(identify){

#if(is.null(BSData$genes)){

#stop("No BSData$genes object found!")

#}

id=identify(0.5*(x+y),y-x, x,n=1, plot=FALSE)

points(0.5*(x[id]+y[id]),y[id]-x[id], pch="X", col="red")

info=getProbeInfo(probes[id])

if(nrow(info)!=0){

print(info)

text(20,2,paste("Probe ID: ", probes[id],sep=""))
text(20,1.8,paste("Gene ID: ", as.character(info[,2]), sep=""))
text(20,1.6,paste("Name:", as.character(info[,1]),sep=""))
text(20,1.4,paste("Value from SAM 1: ", round(x[id],2), sep=""))
text(20,1.2,paste("Value from SAM 2: ", round(y[id],2), sep=""))
}

else{
info=spottypes[spottypes$ID==probes[id],]

print(info)

text(20,2,paste("Probe ID: ", probes[id],sep=""))
text(20,1.8,paste("Control Type: ", as.character(info[,1]),sep=""))
text(20,1.4,paste("Value from SAM 1: ", round(x[id],2), sep=""))
text(20,1.2,paste("Value from SAM 2: ", round(y[id],2), sep=""))

}


}

}




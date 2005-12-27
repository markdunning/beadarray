"plotMA.beads" <-
function(BSData, array1, array2=0, identify=FALSE, label=FALSE, highlight=NULL, log=TRUE, main=NULL){

  if(class(BSData) == "BeadSummaryList"){
    stop("BeadSummaryList object required!")
  }

probes = sort(unique(BSData$probeID[BSData$probeID[,1] > 0,1]))

if(array2!=0){

if(log){

x = log2(BSData$R[,array1])
y = log2(BSData$R[,array2])

}

else{

x = BSData$R[,array1]
y = BSData$R[,array2]

}


}

else{
x = log2(BSData$R[,array1])
y = log2(BSData$G[,array1])

}




plot(0.5*(x+y),y-x,pch=16,cex=0.4,col="black", ylim=range(-2,2),xlim=range(5,25), main=main) 


abline(h=c(-1,0,1),lty=c(2,1,2)) 
  

if(label){

if(is.null(spottypes)){

stop("spottypes object not found")

}


for(i in 1:nrow(spottypes)){

x1 = x[BSData$probeID[,1]==spottypes[i,2]]
y1 = y[BSData$probeID[,1]==spottypes[i,2]]

points(0.5*(x1+y1),y1-x1, col=as.character(spottypes[i,4]))

types = as.character(unique(spottypes[,1]))

cols = as.character(unique(spottypes[,4]))

legend(5,2, types,cols,cex=0.8)

}

}

if(length(highlight!=0)){

for(i in 1:length(highlight)){

x1 = x[BSData$probeID[,1]==highlight[i]]
y1 = y[BSData$probeID[,1]==highlight[i]]


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


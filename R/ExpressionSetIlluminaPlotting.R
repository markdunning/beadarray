plotMAXY <- function(exprs, vec, log = TRUE, labels=vec,label=FALSE,ma.ylim=2,sampleSize=NULL,...){


  mat <- matrix(c(0,1,0,0.04, 0,1,0.96,1, 0,0.04,0.04,0.96,
                0.96,1,0.04,0.96, 0.04,0.96,0.04,0.96), byrow = T, ncol= 4)

close.screen(all=TRUE)

  split.screen(mat)

  split.screen(figs = c(length(vec), length(vec)), screen = 5)
  
  for(i in 1:length(vec)){
    for(j in 1:length(vec)){
      screen(((i-1)*length(vec))+j+5)

      par(mar = c(0.3,0.3,0.3,0.3), cex.axis = 0.7)
      if(i == j){
#        plot(0, col.axis = "white", cex = 0, col.lab = "white", tcl = -0, xlab = "", ylab = "")
        plot(0, axes = TRUE, type = "n", tcl = -0, col.axis = "white")
        text(1.0,0, labels = labels[i], cex=1)
      }
      else if(j < i){
        plotXY.beads(exprs, array1 = vec[i], array2 = vec[j], log = log, xaxt = "n", yaxt = "n", label=label,sampleSize=sampleSize)
        if(i == length(vec)){
          axis(1)
          }
        if(j == 1){
          axis(2)
        }
      }
      else{
        plotMA.beads(exprs, array1 = vec[i], array2 = vec[j], log = log, xaxt = "n", yaxt = "n", label=label,ma.ylim=ma.ylim, sampleSize=sampleSize)
        if(i == 1){
          axis(3)
        }
        if(j == length(vec)){
          axis(4)   
        }
      }
    }
  }
}


"plotMA.beads" <-
function(exprs, array1, array2=0, identify=FALSE, label=FALSE, highlight=NULL, log=TRUE, main=NULL, ma.ylim=2,sampleSize=NULL,...){
exprs=as.matrix(exprs)

if(array2!=0){

if(log){

x = 0.5*(log2(exprs[,array1]) + log2(exprs[,array2]))
y = log2(exprs[,array1])- log2(exprs[,array2])
  

}

else{

x = 0.5*exprs[,array1] + exprs[,array2]
y = exprs[,array1]- exprs[,array2]
}


}

else{
x = log2(exprs[,array1])
y = log2(exprs[,array1])

}

if(!is.null(sampleSize)){
 s = sample(1:length(x), sampleSize)
 x=x[s]
 y=y[s]
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

"plotOnSAM" <-
function(values, mx=max(values, na.rm=TRUE), scale=max(values, na.rm=TRUE),min=0, main=NULL, label=TRUE, missing_arrays=NULL, colour=TRUE){

len = 96

xmax = 1565
ymax = 1565

if(!is.null(missing_arrays)){

values = vector(length=96)

i = 1:96

values[i[-missing_arrays]]=v

values[missing_arrays] = NA

}

par(mfrow=c(1,2))


plot(1:len, values,  ylim=range(min,mx),type="l",  xlab="Array Index", main=main)

abline(v=c(12,24,36,48,60,72,84), lty=3)

if(label){

text(1:len, values, label=1:len)

}


plot(1:xmax*12, xlim=range(0:xmax*12),ylim=range(0:ymax*8), type="n", xlab=" ", ylab=" " ,main=main,xaxt="n", yaxt="n")

ys = c(ymax/2, 0, 0, ymax/2, ymax, ymax)

for(i in 1:8){

xs = c(0,xmax/4, 0.75*xmax, xmax, 0.75*xmax, xmax/4)

for(j in 1:12){

array_index =  (12 * (8-i)) + j

control_intensity = values[array_index]

if(is.na (control_intensity)){
polygon(xs,ys, col=rgb(1,1,1))
}
else{

if(colour){ polygon(xs,ys, col=rgb(control_intensity/scale,0,0))}
else{ polygon(xs,ys, col=gray(control_intensity/scale))}


}

if(label){
text(xs[3], ys[4], array_index)
}

#text(xs[3], ys[4], BLData$other$SAMPLE[1,array_index])

xs = xs + xmax + 20

}

ys = ys + ymax + 30


}



}




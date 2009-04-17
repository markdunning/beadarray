plotMAXY <- function(exprs, arrays, log = TRUE, genesToLabel=NULL,labels=colnames(exprs)[arrays],labelCol="red", labelpch=16,foldLine=2,sampleSize=NULL,...){


  mat <- matrix(c(0,1,0,0.04, 0,1,0.96,1, 0,0.04,0.04,0.96,
                0.96,1,0.04,0.96, 0.04,0.96,0.04,0.96), byrow = T, ncol= 4)

close.screen(all=TRUE)

  split.screen(mat)

  split.screen(figs = c(length(arrays), length(arrays)), screen = 5)
  
  for(i in 1:length(arrays)){
    for(j in 1:length(arrays)){
      screen(((i-1)*length(arrays))+j+5)

      par(mar = c(0.3,0.3,0.3,0.3), cex.axis = 0.7)
      if(i == j){
#        plot(0, col.axis = "white", cex = 0, col.lab = "white", tcl = -0, xlab = "", ylab = "")
        plot(0, axes = TRUE, type = "n", tcl = -0, col.axis = "white")
        text(1.0,0, labels = labels[i], cex=1)
      }
      else if(j < i){
        plotXY(exprs, array1 = arrays[i], array2 = arrays[j], xaxt = "n", yaxt = "n", log=log,genesToLabel=genesToLabel,foldLine=foldLine, labelCol=labelCol,labelpch=labelpch,sampleSize=sampleSize)
        if(i == length(arrays)){
          axis(1)
          }
        if(j == 1){
          axis(2)
        }
      }
      else{
        plotMA(exprs, array1 = arrays[i], array2 = arrays[j], xaxt = "n", yaxt = "n", log=log,genesToLabel=genesToLabel,foldLine=foldLine, labelCol=labelCol,labelpch=labelpch,sampleSize=sampleSize)
        if(i == 1){
          axis(3)
        }
        if(j == length(arrays)){
          axis(4)   
        }
      }
    }
  }
}


"plotMA" <- function(exprs, array1=1, array2=2, genesToLabel=NULL, labelCol="red", foldLine=2, log=TRUE,labelpch=16,ma.ylim=2,sampleSize=NULL,...){
    
    #try to catch any Limma objects and tell the user.
    if(class(exprs)[1] %in% c("RGList", "MAList")) 
        stop("\nIt appears you are trying to use the plotMA() function on a Limma object, but plotMA() is currently masked by beadarray\n\nIf you wish to use the Limma function, you can either call it directly using:\n\t\"limma::plotMA()\"\nor detach the beadarray package using:\n\t\"detach(package:beadarray)\"\n")
    
    exprs=as.matrix(exprs)

    if(log) exprs=log2(exprs)

    if(array2!=0 && array2!=array1){
        x = 0.5*(exprs[,array1] + exprs[,array2])
        y = exprs[,array1]- exprs[,array2]
    }
    else{
        #x = log2(exprs[,array1])
        #y = log2(exprs[,array1])
        stop("\'array1\' and \'array2\' must be different")
    }

    if(!is.null(sampleSize)){
        s = sample(1:length(x), sampleSize)
        x=x[s]
        y=y[s]
    }	
                                                                                                                                       
    naInd = intersect(which(!is.na(x)),which(!is.na(y)))
    naInd = intersect(naInd, which(!(is.infinite(x))))
    naInd = intersect(naInd, which(!(is.infinite(y))))
 	
    smoothScatter(x[naInd],y[naInd], pch=16,cex=0.4, ylim=range(ma.ylim,-ma.ylim), xlab = "", ylab = "", ...) 

    abline(h=c(-log2(foldLine),0,log2(foldLine)),lty=c(2,1,2)) 
  
    if(!is.null(genesToLabel)){
        index = which(rownames(exprs) %in% genesToLabel)
        points(x[index], y[index], col=labelCol, pch=labelpch)
    }
}



"plotXY" <-
function(exprs,array1=1, array2=2, genesToLabel=NULL, labelCol="red", log=TRUE,labelpch=16,foldLine=2,sampleSize=NULL,...){

#XY plot of either two samples against each other, or red and green channels of one channel

if(log) exprs=log2(exprs)
  
exprs=as.matrix(exprs)

  if (array2!=0 && array2!=array1){

 
      x = exprs[,array1]
      y = exprs[,array2]


      xmax = 16
      xbox=18
      yspacing=0.3

                 

   
  }
  else{
#      x = log2(exprs[,array1])
#      y = log2(exprs[,array1])
#
#      xmax = 16
#      xbox=18
#      yspacing=0.3
      stop("\'array1\' and \'array2\' must be different")
    }

  
 if(!is.null(sampleSize)){
 s = sample(1:length(x), sampleSize)
 x=x[s]
 y=y[s]
 }
	

 
    naInd = intersect(which(!is.na(x)),which(!is.na(y)))
  naInd = intersect(naInd, which(!(is.infinite(x))))
  naInd = intersect(naInd, which(!(is.infinite(y))))
                                                                                                                                       
                                                                                                                         

  smoothScatter(x[naInd],y[naInd], xlim=range((max(0,min(x),na.rm=TRUE)),16), xlab = "", ylab = "", pch = 16, cex = 0.4, ...)
  abline(log2(foldLine), 1, lty=2)
  abline(-log2(foldLine),1,lty=2)
  
  if(!is.null(genesToLabel)){

    index = which(rownames(exprs) %in% genesToLabel)

    points(x[index], y[index], col=labelCol, pch=labelpch)
  }

}

"plotOnSAM" <-
function(values, mx=max(values, na.rm=TRUE), scale=max(values, na.rm=TRUE),min=0, main=NULL, label=TRUE, missing_arrays=NULL, colour=TRUE){

len = 96

xmax = 1565
ymax = 1565

if(!is.null(missing_arrays)){

values = vector(length=96)

#i = 1:96

#values[i[-missing_arrays]]=v

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




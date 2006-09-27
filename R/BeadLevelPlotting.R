"plotBeadLocations" <- function(BLData, ProbeIDs=NULL, BeadIDs=NULL, array, SAM=FALSE,...){

xmax = max(BLData@GrnX[,array])
ymax = max(BLData@GrnY[,array])

plot(1, xlim=range(0:xmax), ylim=range(0:ymax) ,type="n", new=TRUE,...)

if(!(is.null(ProbeIDs))){

xs = BLData@GrnX[which(BLData@ProbeID[,array] %in% ProbeIDs),array]
ys = BLData@GrnY[which(BLData@ProbeID[,array] %in% ProbeIDs),array]

}
else{

xs = BLData@GrnX[BeadIDs,array]
ys = BLData@GrnY[BeadIDs,array]
 


}

if(SAM) {

yy = c(ymax/2, 0, 0, ymax/2, ymax, ymax)

xx = c(0,xmax/4, 0.75*xmax, xmax, 0.75*xmax, xmax/4)

polygon(xx, yy)

}

else polygon(x=c(0,0,xmax,xmax), y=c(0, ymax, ymax,0))


points(xs, ys,...)


}



plotBeadIntensities=function(BLData, ProbeIDs, arrays, log=FALSE, n=3, ProbeCols=NULL,ylim=NULL,...){



nplots  = length(ProbeIDs)*length(arrays)

if(is.null(ProbeCols)){

ProbeCols = rainbow(n=length(ProbeIDs), start=0, end=5/6)

}

if(is.null(ylim)) ylim=range(5,12)

plot(4, xlim=range(0,nplots), ylim=ylim, type="n", axes=FALSE, xlab="", ylab="log2 intensities")
count=1


for(i in 1:length(arrays)){
j=1
for(j in 1:length(ProbeIDs)){

I = getProbeIntensities(BLData, ProbeIDs=ProbeIDs[j], array=arrays[i], log=TRUE)

o = findBeadStatus(BLData, probes = ProbeIDs[j], array=arrays[i], log=TRUE, n=n)

o = log2(BLData@G[o,arrays[i]])

boxplot(I, at=count-0.5, add=TRUE, axes=FALSE, col=ProbeCols[j])

points(x=rep(count-0.5, length(o)), y = o, pch=16, col="red")

count = count+1

}

if(i!=length(arrays)) abline(v=count-1)

}
axis(2)
box()

}




"imageplot"<-function(BLData, array = 1, nrow = 18, ncol = 2,
                        low = NULL, high = NULL, ncolors = 123, whatToPlot ="G",zlim=NULL,...){

  par(mar = c(2,1,1,1), xaxs = "i")
  
#Not needed since the co-ords are automatically scaled to zero now  
#  xs <- floor(BLData@GrnX[,array] - min(BLData@GrnX[,array]))
#  ys <- floor(BLData@GrnY[,array] - min(BLData@GrnY[,array]))


  data = slot(BLData, whatToPlot)

  if (is.character(low)) 
    low <- col2rgb(low)/255
  if (is.character(high)) 
    high <- col2rgb(high)/255
  if (!is.null(low) && is.null(high)) 
    high <- c(1, 1, 1) - low
  if (is.null(low) && !is.null(high)) 
    low <- c(1, 1, 1) - high

  if (is.null(low)) 
    low <- c(1, 1, 1)
  if (is.null(high)) 
    high <- c(0, 1, 0)
  
  col <- rgb(seq(low[1], high[1], len = ncolors), seq(low[2], 
        high[2], len = ncolors), seq(low[3], high[3], len = ncolors))

  xs <- floor(BLData@GrnX[,array])
  ys <- floor(BLData@GrnY[,array])

  xgrid <- floor(seq(0, max(xs), by = max(xs)/ncol))
  ygrid <- floor(seq(0, max(ys), by = max(ys)/nrow))

  imageMatrix <- matrix(ncol = ncol, nrow = nrow)

  for(i in 1:ncol){
    idx = which((xs > xgrid[i]) & (xs < xgrid[i+1]))
    fground = data[idx,array]
    yvalues = ys[idx]
#    yvalues = BLData@GrnY[idx,array]

    out <- .C("BLImagePlot", length(fground), as.double(fground), as.double(yvalues), as.integer(ygrid),
              result = double(length = nrow), as.integer(nrow), PACKAGE = "beadarray")

    imageMatrix[,i] <- rev(out$result)
  }

  imageMatrix = t((imageMatrix))
  image(x = c(0:ncol), z = imageMatrix,  xaxt = "n", yaxt = "n", col = col,...)
}



"SAMSummary" <-
function(BLData, mode="outliers", missing_arrays=NULL, colour=TRUE, scale = NULL,...){


#mode=screenSetup(BLData)
  
#Setup up outliers

split.screen(c(1,2))
screen(1)

if(mode=="outliers") title="Number of Outliers"
if(mode=="fg") title="Median Foreground"
if(mode=="bg") title="Median Background"



values = vector(length=96)

if(mode == "outliers"){

  o = list(length=96)

  
  for(i in 1:96){

    o[[i]] = findAllOutliers(BLData, array=i)

    values[i] = length(o[[i]])

  }

}
  
if(mode == "fg"){

 for(i in 1:96){

   values[i] = median(log2(BLData@G[,i]),na.rm=TRUE)
 }
 

}

if(mode == "bg"){

  for(i in 1:96){
    values[i] = median(log2(BLData@Gb[,i]),na.rm=TRUE)
  }
  

}

  
len = 96


if(!is.null(missing_arrays)){

values = vector(length=96)

i = 1:96

values[i[-missing_arrays]]=v

values[missing_arrays] = NA

}


doLoop = TRUE

while(doLoop){

  x=1/13
y=1/9

ys = c(y/2, 0, 0, y/2, y, y)

for(i in 1:8){

xs = c(0,x/4, 0.75*x, x, 0.75*x, x/4)

if(is.null(scale)) scale=max(values)

for(j in 1:12){

array_index =  (12 * (8-i)) + j

control_intensity = values[array_index]

if(is.na (control_intensity)){
polygon(xs,ys, col=rgb(1,1,1))
}
else{

if(colour){ polygon(xs,ys, col=rgb(control_intensity/scale,0,0))}
else{ polygon(xs,ys,col=gray(control_intensity/scale))}

}


#text(xs[3], ys[4], BLData@other@SAMPLE[1,array_index])

xs = xs + 1/12

}

ys = ys + 1/8


}


l=locator(n=1)

x.clicked = strtrim(as.character(l$x),5)
print(l)



y.clicked = strtrim(as.character(l$y),5)

y.clicked = 8-(as.double(y.clicked)*8)

x.clicked = as.double(x.clicked)*12

print(x.clicked)
print(y.clicked)

ArrayClickedOn = ceiling(y.clicked)*12 + ceiling(x.clicked) - 12

print(ArrayClickedOn)

screen(2, new=TRUE)

plot(1:10, type="n", axes=FALSE, xlab="", ylab="")

if(mode == "outliers"){


os= o[[ArrayClickedOn]]

plotBeadLocations(BLData, BeadIDs=os, array=ArrayClickedOn, new=TRUE, main=as.character(BLData@targets[ArrayClickedOn,2]))

screen(1)
}

if(mode == "fg"){

  imageplot(BLData, array=ArrayClickedOn, whatToPlot="G", nrow=50, ncol=50, high="red", low="yellow",main=as.character(BLData@targets[ArrayClickedOn,2]),...)

  screen(1)
}

if(mode == "bg"){

  imageplot(BLData, array=ArrayClickedOn, whatToPlot="Gb", nrow=50, ncol=50, high="red", low="yellow",main=as.character(BLData@targets[ArrayClickedOn,2]),...)
  screen(1)
  
}

}

}

"BeadChipSummary" <-
function(BLData, mode="outliers",colour=TRUE, scale = NULL,...){


  
split.screen(c(1,2))
screen(1)

if(mode=="outliers") title="Number of Outliers"
if(mode=="fg") title="Median Foreground"
if(mode=="bg") title="Median Background"

  
values = vector(length=12)

if(mode == "outliers"){

  o = list(length=12)

  
  for(i in 1:12){

    o[[i]] = findAllOutliers(BLData, array=i)

    values[i] = length(o[[i]])

  }

}
  
if(mode == "fg"){

  values = apply(log2(BLData@G), 2, median)

}


if(mode == "bg"){

  values = apply(log2(BLData@G), 2, median)

}

  

len = 12




doLoop = TRUE

while(doLoop){

  y=1
ys = c(y, y, y-y/12, y - y/12)

  

for(i in 1:12){
control_intensity = values[i]


xs=c(0.1,1,1,0.1)
if(is.null(scale)) scale= max(values)
polygon(xs,ys, col=rgb(control_intensity/scale,0,0))


#text(xs[3], ys[4], BLData@other@SAMPLE[1,array_index])

ys = ys - 1/12

}


l=locator(n=1)


y.clicked = strtrim(as.character(l$y),5)

y.clicked = as.double(y.clicked)*12


ArrayClickedOn = 12 - ceiling(y.clicked) + 1

print(ArrayClickedOn)

screen(2, new=TRUE)

plot(1:10, type="n", axes=FALSE, xlab="", ylab="")

if(mode == "outliers"){


os= o[[ArrayClickedOn]]

plotBeadLocations(BLData, BeadIDs=os, array=ArrayClickedOn, new=TRUE, main=as.character(BLData@targets[ArrayClickedOn,2]))

screen(1)
}

if(mode == "fg"){

  imageplot(BLData, array=ArrayClickedOn, whatToPlot="G", nrow=100, ncol=50, high="red", low="yellow",main=as.character(BLData@targets[ArrayClickedOn,2]))

  screen(1)
}

if(mode == "bg"){

  imageplot(BLData, array=ArrayClickedOn, whatToPlot="Gb", nrow=100, ncol=50, high="red", low="yellow",main=as.character(BLData@targets[ArrayClickedOn,2]))
  screen(1)
  
}

}

}
 






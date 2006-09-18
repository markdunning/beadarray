"plotBeadLocations" <- function(BLData, ProbeIDs=NULL, beadIDs=NULL, array, label=FALSE,...){

dims = getDimensions(BLData, array)
xmax = dims$xmax
ymax = dims$ymax


if(!is.null(ProbeIDs)){

xs = ys = NULL

for(i in 1:length(ProbeIDs)){

xs = c(xs, getProbeCoords(BLData, ProbeIDs[i], array)[,1])
ys = c(ys, getProbeCoords(BLData, ProbeIDs[i], array)[,2])


}

coords = matrix(nrow=length(xs), ncol=2)

#standardise the coordinates to start from 0

coords[,1] = xs 
coords[,2] = ys 


plotCoords(xmax, ymax, coords, label=label,...)


}

else{

coords=matrix(nrow=length(beadIDs), ncol=2)

coords[,1] = BLData$GrnX[beadIDs, array] - min(BLData$GrnX[,array])
coords[,2] = BLData$GrnY[beadIDs, array] - min(BLData$GrnY[,array])


plotCoords(xmax, ymax, coords,label=label,...)

}

}

plotBeadIntensities=function(BLData, ProbeIDs, arrays, log=FALSE, n=3, ProbeCols=NULL,ylim=NULL,...){



nplots  = length(ProbeIDs)*length(arrays)

if(is.null(ProbeCols)){

ProbeCols = rainbow(n=length(ProbeIDs), start=0, end=5/6)

}

if(is.null(ylim)) ylim=range(5,12)

plot(4, xlim=range(0,nplots), ylim=ylim, type="n")
count=1


for(i in 1:length(arrays)){

for(j in 1:length(ProbeIDs)){

I = getProbeIntensities(BLData, ProbeIDs=ProbeIDs[j], array=arrays[i], log=TRUE)

o = findBeadStatus(BLData, probes = ProbeIDs[j], array=arrays[i], log=TRUE, n=n)

o = log2(BLData$G[o,arrays[i]])

boxplot(I, at=count-0.5, add=TRUE, axes=FALSE, col=ProbeCols[j])

points(x=rep(count-0.5, length(o)), y = o, pch=16, col="red")

count = count+1

}

abline(v=count-1)

}


}

"plotCoords" <-
function(xmax, ymax, coords, label,...){

plotArray(xmax, ymax,...)

#plotArray draws out the hexagonal shape of the array and divides into 8 sections


if(label){



#On actual arrays, the y co-ordinates are arranged with y=0 at the top of the array, 
#whereas in R we plot with y=0 at the bottom of the plot. Therefore we transform the 
#y co-ordinates before plotting

text(coords[,1], (ymax - coords[,2]),labels = seq(1:length(coords[,1])),...)

}

else{

points(coords[,1], (ymax - coords[,2]),...)

}


}

"plotArray" <-
function(xmax, ymax, SAM=TRUE,...){

  plot(0,type="n",xlab="",ylab="",ylim=range(0,ymax), xlim=range(0,xmax))

if(!SAM){

  polygon(x=c(0,0,xmax,xmax), y=c(0, ymax, ymax,0))

}


else{





#Plots out the hexagonal arrangement of an array array 

ys = c(ymax/2, 0, 0, ymax/2, ymax, ymax)

xs = c(0,xmax/4, 0.75*xmax, xmax, 0.75*xmax, xmax/4)

polygon(xs, ys)

#Plot lines through centre of hexagon

abline(v=xmax/2)

abline(h=ymax/2)


#Plot circle of radius xmax/4 at the centre
r= xmax/4

theta<-seq(0, 2*pi, length=100)

lines(r*cos(theta)+xmax/2, r*sin(theta)+ymax/2)

}

}


"imageplot"<-function(BLData, array = 1, nrow = 18, ncol = 2,
                        low = NULL, high = NULL, ncolors = 123, whatToPlot ="G"){

  par(mar = c(2,1,1,1), xaxs = "i")
  
#Not needed since the co-ords are automatically scaled to zero now  
#  xs <- floor(BLData$GrnX[,array] - min(BLData$GrnX[,array]))
#  ys <- floor(BLData$GrnY[,array] - min(BLData$GrnY[,array]))

  idx = which(names(BLData) == whatToPlot)

  data = BLData[[idx]]		

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

  xs <- floor(BLData$GrnX[,array])
  ys <- floor(BLData$GrnY[,array])

  xgrid <- floor(seq(0, max(xs), by = max(xs)/ncol))
  ygrid <- floor(seq(0, max(ys), by = max(ys)/nrow))

  imageMatrix <- matrix(ncol = ncol, nrow = nrow)

  for(i in 1:ncol){
    idx = which((xs > xgrid[i]) & (xs < xgrid[i+1]))
    fground = data[idx,array]
    yvalues = ys[idx]
#    yvalues = BLData$GrnY[idx,array]

    out <- .C("BLImagePlot", length(fground), as.double(fground), as.double(yvalues), as.integer(ygrid),
              result = double(length = nrow), as.integer(nrow), PACKAGE = "beadarray")

    imageMatrix[,i] <- rev(out$result)
  }

  imageMatrix = t((imageMatrix))
  image(x = c(0:ncol), z = imageMatrix,  xaxt = "n", yaxt = "n", col = col)
}


screenSetup=function(BLData){

  pushViewport(viewport(x=0, width=0.1, height=1,name="mode.selector",just=c("left")))
  grid.rect(y=1, height=0.1)
  grid.text(y=unit(1, "npc")-unit(1, "char"), x=unit(1,"npc")-unit(1, "char"),label="O")
  grid.rect(y=0.9, height=0.1)
  grid.text(y=unit(0.9, "npc")-unit(1, "char"), x=unit(1,"npc")-unit(1, "char"),label="Fg")
  grid.rect(y=0.8, height=0.1)
  grid.text(y=unit(0.8, "npc")-unit(1, "char"), x=unit(1, "npc")-unit(1, "char"),label="Bg")
  
  upViewport()
  pushViewport(viewport(x=0.1, width=1, height=1, name="main",default.units="npc",just=c("left")))
  grid.rect()


  l=grid.locator(unit="npc")

  y = as.double(strtrim(as.character(l$y),5))

  
  if(y > 0.9 & y<1){

    m = "outliers"

  }

  if(y > 0.8 & y<0.9){

    m = "fg"

  }

  if(y>0.7 & y<0.8){

    m="bg"

  }

  
  print(m)
  
  
  m


}


"SAMSummary" <-
function(BLData, mode="outliers"){


#mode=screenSetup(BLData)
  
#Setup up outliers

grid.rect(gp=gpar(fill="white"))
  
pushViewport(viewport(width=0.48, height=0.9,name="array.selector",  just="right"))
grid.rect()

if(mode=="outliers") title="Number of Outliers"
if(mode=="fg") title="Median Foreground"
if(mode=="bg") title="Median Background"

grid.text(y=unit(1, "npc")+unit(1, "char"), label=title)
upViewport()

pushViewport(viewport(width=0.48, height=0.9,name="array.viewer",  just="left"))


upViewport()



seekViewport("array.selector")

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

   values[i] = median(log2(BLData$R[,i]),na.rm=TRUE)
 }
 

}

if(mode == "bg"){

  for(i in 1:96){
    values[i] = median(log2(BLData$Rb[,i]),na.rm=TRUE)
  }
  

}

  
  

len = 96


if(!is.null(missing_arrays)){

values = vector(length=96)

i = 1:96

values[i[-missing_arrays]]=v

values[missing_arrays] = NA

}

x=1/13
y=1/9

ys = c(y/2, 0, 0, y/2, y, y)

for(i in 1:8){

xs = c(0,x/4, 0.75*x, x, 0.75*x, x/4)

for(j in 1:12){

array_index =  (12 * (8-i)) + j

control_intensity = values[array_index]

if(is.na (control_intensity)){
grid.polygon(xs,ys, col=rgb(1,1,1))
}
else{

if(colour){ grid.polygon(xs,ys, gp=gpar(fill=rgb(control_intensity/scale,0,0)),default.units="npc")}
else{ grid.polygon(xs,ys, gp=gpar(col=gray(control_intensity/scale)))}


}

if(label){
text(xs[3], ys[4], array_index)
}

#text(xs[3], ys[4], BLData$other$SAMPLE[1,array_index])

xs = xs + 1/12

}

ys = ys + 1/8


}

doLoop = TRUE

while(doLoop){

l=grid.locator(unit="npc")

x.clicked = strtrim(as.character(l$x),5)
print(l)

y.clicked = strtrim(as.character(l$y),5)

y.clicked = 8-(as.double(y.clicked)*8)

x.clicked = as.double(x.clicked)*12

print(x.clicked)
print(y.clicked)

ArrayClickedOn = ceiling(y.clicked)*12 + ceiling(x.clicked) - 12

print(ArrayClickedOn)

seekViewport("array.viewer")

grid.rect(gp=gpar(fill="white"))

grid.text(x=0.5, y=unit(1, "npc")+unit(1, "char"), label="Outlier Locations")

if(mode == "outliers"){


os= o[[ArrayClickedOn]]

#Transform x and y coordinates to range 0,1

xs.npc = (BLData$x[os, ArrayClickedOn] - min(BLData$x[,ArrayClickedOn])) / max(BLData$x[os, ArrayClickedOn])*0.9+0.01

ys.npc = (BLData$y[os, ArrayClickedOn] - min(BLData$y[,ArrayClickedOn])) / max(BLData$y[os, ArrayClickedOn])

grid.points(x=xs.npc, y=ys.npc, pch=16, size=unit(1,"mm"),gp=gpar(col="red"))

seekViewport("array.selector")

}




}

}



"BeadChipSummary" <-
function(BLData, mode="outliers"){

#mode=screenSetup(BLData)
  

pushViewport(viewport(x=0,width=0.3, height=0.90,name="array.selector",default.units="npc",just="left"))
grid.rect()
upViewport()

pushViewport(viewport(x=0.3,width=0.9, height=0.90,name="array.viewer",default.units="npc",just="left"))
grid.rect()
upViewport()



seekViewport("array.selector")

values = vector(length=12)

if(mode == "outliers"){

  o = list(length=12)

  
  for(i in 1:12){

    o[[i]] = findAllOutliers(BLData, array=i)

    values[i] = length(o[[i]])

  }

}
  
if(mode == "fg"){

  values = apply(log2(BLData$R), 2, median)

}

  
  

len = 12


if(!is.null(missing_arrays)){

values = vector(length=96)

i = 1:12

values[i[-missing_arrays]]=v

values[missing_arrays] = NA

}


y=1
ys = c(y, y, y-y/12, y - y/12)

for(i in 1:12){
control_intensity = values[i]


xs=c(0.1,1,1,0.1)

if(colour){ grid.polygon(xs,ys, gp=gpar(fill=rgb(control_intensity/scale,0,0)),default.units="npc")}
else{ grid.polygon(xs,ys, gp=gpar(col=gray(control_intensity/scale)))}





#text(xs[3], ys[4], BLData$other$SAMPLE[1,array_index])

ys = ys - 1/12

}



doLoop = TRUE

while(doLoop){

l=grid.locator()


y.clicked = strtrim(as.character(l$y),5)

y.clicked = as.double(y.clicked)*12


ArrayClickedOn = 12 - ceiling(y.clicked) + 1

print(ArrayClickedOn)

seekViewport("array.viewer")

grid.rect(gp=gpar(fill="white"))

grid.text(x=0.5, y=unit(1, "npc")+unit(1, "char"), label="Outlier Locations")

if(mode == "outliers"){


os= o[[ArrayClickedOn]]

#Transform x and y coordinates to range 0,1

xs.npc = (BLData$x[os, ArrayClickedOn] - min(BLData$x[,ArrayClickedOn])) / max(BLData$x[os, ArrayClickedOn]+0.1)*0.7 + 0.01

ys.npc = (BLData$y[os, ArrayClickedOn] - min(BLData$y[,ArrayClickedOn])) / max(BLData$y[os, ArrayClickedOn])

grid.points(x=xs.npc, y=ys.npc, pch=16, size=unit(0.5,"mm"),gp=gpar(col="red"))

seekViewport("array.selector")

}




}

}












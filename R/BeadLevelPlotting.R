boxplotBeads = function(BLData, whatToPlot="G", arrays=NULL,
                                  log=TRUE, n=3, varwidth=TRUE,  ...) {
  tmp = list()
  arraynms = arrayNames(BLData)
  narrays = length(arraynms)
  if(is.null(arrays))  # plot all arrays
    arrays = 1:narrays
  for(i in arrays)
      tmp[[arraynms[i]]] = getArrayData(BLData,array=i, which=whatToPlot, log=log, n=n)
  boxplot(tmp,varwidth=varwidth,...)
}

plotRG = function(BLData, ProbeIDs=NULL, BeadIDs=NULL, log=TRUE, arrays=1,
                   xlim=c(8,16), ylim=c(8,16), xlab="G intensities",
                   ylab="R intensities", main=arrayNames(BLData)[arrays],
                   smooth=TRUE, cols=NULL, ...) {
  arraynms = arrayNames(BLData)
  narrays = length(arrays)
  if(length(arrays)==1 & is.null(ProbeIDs) & is.null(BeadIDs)) {
     if(smooth)
       smoothScatter(getArrayData(BLData, which="G", array=arrays, log=log),
            getArrayData(BLData, which="R", array=arrays, log=log),
            xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab,
            main=main, ...)
     else
       plot(getArrayData(BLData, which="G", array=arrays, log=log),
            getArrayData(BLData, which="R", array=arrays, log=log),
            xlim=xlim, ylim=ylim, new=TRUE, xlab=xlab, ylab=ylab,
            main=main, ...)
  }
  else{
    if(is.null(cols))
      cols = rainbow(narrays)
    plot(1, xlim=xlim, ylim=ylim, type="n", new=TRUE, xlab=xlab, ylab=ylab, main=main, ...)
     for(i in arrays) {
           if(!is.null(ProbeIDs))
             for(j in 1:length(ProbeIDs)) {
               points(getArrayData(BLData, which="G", array=i, log=log)[which(BLData[[i]]$ProbeID %in% ProbeIDs[j])],
                      getArrayData(BLData, which="R", array=i, log=log)[which(BLData[[i]]$ProbeID %in% ProbeIDs[j])], col = cols[which(arrays %in% i)], ...)
                }
            if (!is.null(BeadIDs)) 
                points(getArrayData(BLData, which = "G", array = i, 
                  log = log)[BeadIDs], getArrayData(BLData, which = "R", 
                  array = i, log = log)[BeadIDs], col = cols[which(arrays %in% i)],
                  ...)
        }
    }
}

plotBeadLocations = function(BLData, ProbeIDs=NULL, BeadIDs=NULL, array=1, SAM=FALSE, xlab="x-coordinate", ylab="y-coordinate", main=paste("Bead", ProbeIDs, "locations"),...){

xmax = max(BLData[[array]]$GrnX)
ymax = max(BLData[[array]]$GrnY)

plot(1, xlim=range(0:xmax), ylim=range(0:ymax) ,type="n", new=TRUE, xlab=xlab, ylab=ylab,main=main,...)

if(!(is.null(ProbeIDs))){

xs = BLData[[array]]$GrnX[which(BLData[[array]]$ProbeID %in% ProbeIDs)]
ys = ymax-BLData[[array]]$GrnY[which(BLData[[array]]$ProbeID %in% ProbeIDs)]

}
else{

xs = BLData[[array]]$GrnX[BeadIDs]
ys = ymax-BLData[[array]]$GrnY[BeadIDs]
 
}

if(SAM) {

yy = c(ymax/2, 0, 0, ymax/2, ymax, ymax)

xx = c(0,xmax/4, 0.75*xmax, xmax, 0.75*xmax, xmax/4)

polygon(xx, yy)

}

else polygon(x=c(0,0,xmax,xmax), y=c(0, ymax, ymax,0))


points(xs, ys,...)


}


plotBeadIntensities = function(BLData, ProbeIDs, arrays, log=FALSE, whatToPlot="G", ProbeCols=NULL, ylim=NULL,...){

nplots  = length(ProbeIDs)*length(arrays)

if(is.null(ProbeCols)){

ProbeCols = rainbow(n=length(ProbeIDs), start=0, end=5/6)

}

if(is.null(ylim)) {
  if(log==TRUE)
    ylim=c(0,16)
  else
    ylim=c(1,2^16)
}

plot(4, xlim=range(0,nplots), ylim=ylim, type="n", axes=FALSE, xlab="", ylab="intensities",...)
count=1


for(i in 1:length(arrays)){
for(j in 1:length(ProbeIDs)){
I = getProbeIntensities(BLData, ProbeIDs=ProbeIDs[j], array=arrays[i], log=log, which=whatToPlot)
  if (length(I) > 1) {
      sel = is.finite(I) & !is.na(I)
      boxplot(I[sel], at = count - 0.5, add = TRUE, axes = FALSE, col = ProbeCols[j])
  }
  else {
    cat("\nNo bead with ID", ProbeIDs[j], "on array", i, "\n")
  }
#points(x=rep(count-0.5, length(o)), y = o, pch=16, col=ProbeCols[j])

count = count+1

}

if(i!=length(arrays)) abline(v=count-1)

}
axis(2)
box()

}




imageplot = function(BLData, array = 1, nrow = 100, ncol = 100,
                     low = NULL, high = NULL, ncolors = 123,
                     whatToPlot ="G", log=TRUE, n=3, zlim=NULL,
                     main=whatToPlot,...){

  par(mar = c(2,1,1,1), xaxs = "i")
  
#Not needed since the co-ords are automatically scaled to zero now  
#  xs = floor(BLData@GrnX[,array] - min(BLData@GrnX[,array]))
#  ys = floor(BLData@GrnY[,array] - min(BLData@GrnY[,array]))

  whatToPlot = match.arg(whatToPlot, choices=c("G", "Gb", "R", "Rb", "wtsG", "wtsR", "residG", "residR", "M", "residM", "A"))
  if((whatToPlot=="R" | whatToPlot=="residR" | whatToPlot=="M" | whatToPlot=="residM" | whatToPlot=="A") & BLData@arrayInfo$channels!="two")
    stop(paste("Need two-channel data to plot", whatToPlot, "values"))
                                          
  data = getArrayData(BLData, which=whatToPlot, array=array, log=log) 
  ind = is.na(data) | is.infinite(data)
  if(sum(ind)>0) {
    cat(paste("Warning:", sum(ind), "NA, NaN or Inf values, which will be set to zero.\nCheck your data or try setting log=\"FALSE\"\n"))
    data[ind] = 0
    rm(ind)
  }
  if (is.character(low)) 
    low = col2rgb(low)/255
  if (is.character(high)) 
    high = col2rgb(high)/255
  if (!is.null(low) && is.null(high)) 
    high = c(1, 1, 1) - low
  if (is.null(low) && !is.null(high)) 
    low = c(1, 1, 1) - high

  if (is.null(low)) 
    low = c(1, 1, 1)
  if (is.null(high)) 
    high = c(0, 1, 0)
#  if(whatToPlot=="G" | whatToPlot=="Gb")
    col = rgb(seq(low[1], high[1], len = ncolors), seq(low[2], 
          high[2], len = ncolors), seq(low[3], high[3], len = ncolors))

#  else  # plot in Red colour scheme
#    col = rgb(seq(low[2], high[2], len = ncolors), seq(low[1], 
#          high[1], len = ncolors), seq(low[3], high[3], len = ncolors))

  xs = floor(BLData[[array]]$GrnX)
  ys = floor(BLData[[array]]$GrnY)

  xgrid = floor(seq(0, max(xs), by = max(xs)/ncol))
  ygrid = floor(seq(0, max(ys), by = max(ys)/nrow))

  imageMatrix = matrix(ncol = ncol, nrow = nrow)

  for(i in 1:ncol){
    idx = which((xs > xgrid[i]) & (xs < xgrid[i+1]))
    if(length(idx)>0) {
      fground = data[idx]
      yvalues = ys[idx]
#      yvalues = BLData@GrnY[idx,array]

      out = .C("BLImagePlot", length(fground), as.double(fground), as.double(yvalues), as.integer(ygrid),
              result = double(length = nrow), as.integer(nrow), PACKAGE = "beadarray")
      if (!is.null(zlim)) {
         out$result[!is.na(out$result)] = pmax(zlim[1], out$result[!is.na(out$result)], na.rm=TRUE)
         out$result[!is.na(out$result)] = pmin(zlim[2], out$result[!is.na(out$result)], na.rm=TRUE)
      }
      imageMatrix[,i] = rev(out$result)
    }
  }
 
  imageMatrix = t((imageMatrix))
  if(is.null(zlim)) zlim=range(imageMatrix, na.rm=TRUE)
  image(x = c(0:ncol), z = imageMatrix,  xaxt = "n", yaxt = "n", col = col, main=main,zlim=zlim,...)
}



SAMSummary =
function(BLData, mode="outliers", whatToPlot="G", log=TRUE, n=3, missing_arrays=NULL, colour=TRUE, scale = NULL, low="yellow", high="red",...) {

#mode=screenSetup(BLData)

#Setup up outliers

split.screen(c(1,2))
screen(1)

#if(mode=="outliers") title=paste("Number of Outliers:", whatToPlot)
#if(mode=="intensities") title=paste("Median Intensities:", whatToPlot)

arraynms = arrayNames(BLData)
narrays = length(arraynms)
len = min(96, narrays)                                                                                                           
values = vector(length=96)

if(mode == "outliers"){

  o = list(length=len)
  
  for(i in 1:len){

    o[[i]] = findAllOutliers(BLData, array=i, log=log, which=whatToPlot, n=n,...)

    values[i] = length(o[[i]])

  }

}
  
if(mode == "intensities") {

 for(i in 1:len) {

   values[i] =  median(getArrayData(BLData, array=i, log=log, which=whatToPlot), na.rm=TRUE) # median(log2(BLData@G[,i]),na.rm=TRUE)
 }
 
}
#if(mode == "bg"){
#
#  for(i in 1:96){
#    values[i] = median(log2(BLData@Gb[,i]),na.rm=TRUE)
#  }
#  
#
#}
#len = 96

if(!is.null(missing_arrays)){

values = vector(length=len)

i = 1:len

values[i[-missing_arrays]]=v

values[missing_arrays] = NA

}


doLoop = TRUE
counter = 0
while(doLoop){
counter = counter+1
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
if(counter>1) {
  answer = readline("Type Q to quit, or anything else to continue:  ")
  if(answer=="Q")
    break# doLoop=FALSE
}
cat(paste("\nClick on plot device to select an array to show", mode, "from\n\n"))
l=locator(n=1)

x.clicked = strtrim(as.character(l$x),5)
#print(l)

y.clicked = strtrim(as.character(l$y),5)

y.clicked = 8-(as.double(y.clicked)*8)

x.clicked = as.double(x.clicked)*12

#print(x.clicked)
#print(y.clicked)

ArrayClickedOn = ceiling(y.clicked)*12 + ceiling(x.clicked) - 12

#print(ArrayClickedOn)
cat(paste("Displaying", mode, "from array", arraynms[ArrayClickedOn], "\n\n"))

screen(2, new=TRUE)

plot(1:10, type="n", axes=FALSE, xlab="", ylab="")

if(mode == "outliers"){


os= o[[ArrayClickedOn]]

plotBeadLocations(BLData, BeadIDs=os, array=ArrayClickedOn, SAM=TRUE, new=TRUE, main=paste("Outliers:", whatToPlot))
screen(1)
}

if(mode == "intensities") {

  imageplot(BLData, array=ArrayClickedOn, whatToPlot=whatToPlot, nrow=50, ncol=50, high=high, low=low,main=paste("Intensities:", whatToPlot),...)
  screen(1)
}
}
}

"BeadChipSummary" =
function(BLData, mode="outliers", whatToPlot="G", log=TRUE, n=3, colour=TRUE, scale = NULL,low="yellow", high="red",...){

split.screen(c(1,2))
screen(1)

#if(mode=="outliers") title=paste("Number of Outliers:", whatToPlot)
#if(mode=="intensities") title=paste("Median Intensities:", whatToPlot)
arraynms = arrayNames(BLData)
narrays = length(arraynms)
len = min(12, narrays)                                                          
values = vector(length=len)

if(mode == "outliers"){

  o = list(length=len)

  
  for(i in 1:len){

    o[[i]] = findAllOutliers(BLData, array=i, log=log, n=n, which=whatToPlot)
    values[i] = length(o[[i]])

  }

}

if(mode == "intensities"){

 for(i in 1:len){
   values[i] =  median(getArrayData(BLData, which=whatToPlot, log=log, array=i), na.rm=TRUE) # median(log2(BLData@G[,i]),na.rm=TRUE)
 }
 

}
#if(mode == "fg"){

#  values = getArrayData(BLData, which=whatToPlot, array=array, log=log)
# apply(log2(BLData@G), 2, median)

#}


#if(mode == "bg"){

#  values = apply(log2(BLData@G), 2, median)

#}


#len = 12

doLoop = TRUE
counter = 0
while(doLoop){
counter = counter + 1
  y=1
ys = c(y, y, y-y/12, y - y/12)

  

for(i in 1:len){
control_intensity = values[i]


xs=c(0.1,1,1,0.1)
if(is.null(scale)) scale= max(values)
polygon(xs,ys, col=rgb(control_intensity/scale,0,0))


#text(xs[3], ys[4], BLData@other@SAMPLE[1,array_index])

ys = ys - 1/12

}

if(counter>1) {
  answer = readline("Type Q to quit, or anything else to continue:  ")
  if(answer=="Q")
    break# doLoop=FALSE
}
cat(paste("\nClick on plot device to select an array to show", mode, "from\n\n"))
l=locator(n=1)


y.clicked = strtrim(as.character(l$y),5)

y.clicked = as.double(y.clicked)*12


ArrayClickedOn = 12 - ceiling(y.clicked) + 1

cat(paste("Displaying", mode, "from array", arraynms[ArrayClickedOn], "\n\n"))

screen(2, new=TRUE)

plot(1:10, type="n", axes=FALSE, xlab="", ylab="")

if(mode == "outliers"){


os= o[[ArrayClickedOn]]

plotBeadLocations(BLData, BeadIDs=os, array=ArrayClickedOn, new=TRUE, main=paste("Outliers:", whatToPlot))

screen(1)
}

if(mode == "intensities"){

  imageplot(BLData, array=ArrayClickedOn, whatToPlot=whatToPlot, nrow=100, ncol=50, high=high, low=low, main=paste("Intensities:", whatToPlot),...)

  screen(1)
}
#if(mode == "bg"){

#  imageplot(BLData, array=ArrayClickedOn, whatToPlot=whatToPlot, nrow=100, ncol=50, high="red", low="yellow",main=as.character(arrayNames(BLData)[ArrayClickedOn]))
#  screen(1)
  
#}

}

}

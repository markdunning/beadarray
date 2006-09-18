"displayTIFFImage" <-
function(BLData, array, a=680:720, b=680:720,flip=TRUE,showOutliers=TRUE,locateBeads=FALSE, showUnregistered=FALSE, outliers=NULL){

tif = as.character(BLData$targets$Image1[array])

print(tif)

xt = .Call("readTIFF", tif, PACKAGE = "beadarray")

if(showOutliers & is.null(outliers)){
outliers = findAllOutliers(BLData, array)
}

if(showUnregistered){
u = getUnregisteredBeads(BLData, 1)
}
else u=NULL

if(locateBeads){
showTIFF(xt, BLData, BLData$x[,array], BLData$y[,array], a, b, array=array, outliers=outliers,out=showOutliers,locate="beads",flip=flip, unregs=u)
}
else showTIFF(xt, BLData,BLData$x[,array], BLData$y[,array], a, b, array=array, outliers=outliers,out=showOutliers,flip=flip,unregs=u)


}


"showTIFF" <-
function(xt,BLData,xs,ys,a,b, array, zScale=TRUE, 
                         reverseYaxis = TRUE,
                        width=9, height=9,
                       resc=FALSE,showSpots=TRUE, out=FALSE,outliers, unregs=NULL, contrast=FALSE,locate=FALSE,recenter=FALSE,flip=FALSE){



info=list(width=width,height=height,Yaxis=reverseYaxis,a=a,b=b,scale=zScale,resc=resc,out=out,showSpots=showSpots,
contrast=contrast,locate=locate)


if (flip){

b  = sort(max(ys) - b)

}

  if( ! is.matrix(xt) ){
    warning(" class matrix is required in plot.imageMatrix", call.=FALSE)
    return()
  }
  if( ! all(is.numeric(xt)) ){
    cat(" matrix x is automatically coerced to numeric in plot.imageMatrix \n")
    xt <- as.numeric(xt)
  }

par(ask=TRUE)
x<-xt[a,b]
ap=a
if(zScale){
a=min(a):(max(a)+floor((max(a)-min(a))/4))
x<-xt[a,b]
}
lambda<-max(a)-min(a)
mu<-max(b)-min(b)

 #x=Contrast(x,contrast)
   
  nrOfRows <- dim(x)[1]       ## relates to length of y-axis
  nrOfColumns <- dim(x)[2]    ## relates to length of x-axis
 
  ## reverse axis
  if( reverseYaxis ){
    reverse <- nrOfColumns:1
    x <- x[,reverse]
      }
    
   
x11(width=width,height=height)
      ## just trying to avoid trouble with Windows Screens ...
   if( .Platform$OS.type == "windows" ){
      par(mar=c(3, 2.5, 2.5, 2))
   }
    
colourRamp <- rgb(0,seq(0,1,l=256),0)

 if(resc==TRUE){
 image(a, b, log2(x), xlab="x", ylab="y", col=colourRamp, axes=TRUE)
 }
 if(resc==FALSE){
    for(i in 1:length(a)){
 for(j in 1:length(b)){
  if(x[i,j] < 1){
 x[i,j] <-1
}
}
  }



  image(a, b, log2(x), xlab="x", ylab="y", axes=TRUE,col=colourRamp, zlim=c(0,16))

}
  
  
##Shows spots 
  if(showSpots==TRUE){
  FindIlluminaSpotl(xt,xs,ys,a,b)
  }

##Shows outliers
if(out==TRUE){

print("finding the outliers....")

redspots<-0
bluespots<-0
yellowspots<-0
unregspots<-0


d<-(xs[outliers]< max(ap) & xs[outliers] >min(ap) & ys[outliers]< max(b) & ys[outliers]>min(b))





outliers_to_plot=outliers[d]
dyellow=(BLData$G[outliers[d],array]<0)
dpos=(BLData$G[outliers[d],array]>0)


for(j in outliers_to_plot[dyellow]){
points(xs[j],max(b)-ys[j]+min(b),col="yellow",cex=1,pch=19)
yellowspots<-yellowspots+1
}

if (!missing(unregs)){

u<-(xs[unregs]< max(ap) & xs[unregs] >min(ap) & ys[unregs]< max(b) & ys[unregs]>min(b))

unregs_to_plot = unregs[u]

for(j in unregs_to_plot){

points(BLData$x[j,array],max(b)-BLData$y[j,array]+min(b), col="grey", cex=1, pch=19)


unregspots = unregspots + 1

}


}



out_blue_red=outliers_to_plot[dpos]
test=vector(length=length(out_blue_red))
for(j in 1:(length(out_blue_red))){
test[j]=BLData$G[out_blue_red[j],array]-mean(getProbeIntensities(BLData, BLData$ProbeID[out_blue_red[j],array], array, log=FALSE), na.rm=TRUE)
}

dblue=(test>0)
dred=(test<0)

for(j in out_blue_red[dblue]){
points(BLData$x[j,array],max(b)-BLData$y[j,array]+min(b),col="blue",cex=1,pch=19)

bluespots<-bluespots+1
}

for(j in out_blue_red[dred]){

points(BLData$x[j,array],max(b)-BLData$y[j,array]+min(b),col="red",cex=1,pch=19)
redspots<-redspots+1

}

 
##Scale
if(zScale){
     
      rect(max(ap)+0.5,min(b)-0.5,max(a)+0.5,max(b)+0.5,col="white")
  if(contrast){
  
  }
  for(i in 1:16){
  rect(max(ap)+lambda*0.08,min(b)+mu*(0.4+0.03*(i-1)),max(ap)+0.5+lambda*0.12,min(b)+mu*(0.4+0.03*i),col=colourRamp[16*i],border=FALSE)
segments(max(ap)+0.5+lambda*0.12,min(b)+mu*(0.4+0.03*i),max(ap)+0.5+lambda*0.125,min(b)+mu*(0.4+0.03*i))
text(x=max(ap)+0.5+lambda*0.14,y=min(b)+mu*(0.4+0.03*i),labels=i,cex=0.7)
    } } 

 
 
 ##Legend 
 if(zScale){
if(out==TRUE){
legend(max(ap)+lambda*0.1,min(b)+0.1*mu,c(">mean","<mean","<0", "unreg"),col=c("blue","red","yellow","white"),pch=19,cex=0.7)
print(paste("<mean",redspots)) 
print(paste(">mean",bluespots)) 
print(paste("negative intensity",yellowspots)) 
text(max(ap)+lambda*0.10,min(b)+mu*0.23,lwd=2,bluespots,cex=0.8)
text(max(ap)+lambda*0.10,min(b)+mu*0.2,lwd=2,redspots,cex=0.8)
text(max(ap)+lambda*0.10,min(b)+mu*0.17,lwd=2,yellowspots,cex=0.8)
text(max(ap)+lambda*0.10,min(b)+mu*0.14,lwd=2,unregspots,cex=0.8)

points(max(ap)+lambda*0.16,min(b)+mu*0.23,lwd=2,pch=19,cex=1.5,col="blue")
points(max(ap)+lambda*0.16,min(b)+mu*0.2,lwd=2,pch=19,cex=1.5,col="red")
points(max(ap)+lambda*0.16,min(b)+mu*0.17,lwd=2,pch=19,cex=1.5,col="yellow")
points(max(ap)+lambda*0.16,min(b)+mu*0.14,lwd=2,pch=19,cex=1.5,col="grey")

}
}

if(recenter)
{
locate=FALSE
Recenter(xt,xs,ys,info,outliers, recenter)
}

if(locate=="beads"){
id<-identify(xs,max(b)-ys+min(b),BLData$ProbeID[,array],plot=FALSE,pos=FALSE,n=1)
IdentifyBeads(a+0.5,ap+0.5,b,lambda,mu,id, outliers=outliers, xt=xt, array=array,xs=xs,ys=ys)

}
}
  if(locate=="intensity"){
  loc<-locator(1)
IdentifyIntensity(a+0.5,ap+0.5,b,lambda,mu,loc)
  }
  
  
  
  
}# end of the function definition plot.imageMatrix4

"IdentifyBeads" <-
function(a,ap,b,lambda,mu,id,xt,array,outliers, xs, ys,...){

rect(max(ap)+0.01*lambda,min(b)+0.3*mu,max(a)-0.03*lambda,max(b),col="white",border=FALSE)

if(is.na(match(id,outliers))==FALSE){
outliers[id]<-id
text(max(ap)+lambda*0.16,max(b)-mu*0.05,"Outlier")
text(max(ap)+lambda*0.16,max(b)-mu*0.1,lwd=2,xs[id],cex=0.8)
text(max(ap)+lambda*0.16,max(b)-mu*0.14,lwd=2,ys[id],cex=0.8)
text(max(ap)+lambda*0.16,max(b)-mu*0.18,lwd=2,round(BLData$G[outliers[id],array],2),cex=0.8)
text(max(ap)+lambda*0.16,max(b)-mu*0.22,lwd=2,round(BLData$Gb[outliers[id],array],2),cex=0.8)
text(max(ap)+lambda*0.16,max(b)-mu*0.26,lwd=2,round(mean(getProbeIntensities(BLData,BLData$ProbeID[outliers[id],array],array, log=FALSE)),2),cex=0.8)
text(max(ap)+lambda*0.16,max(b)-mu*0.30,lwd=2,BLData$ProbeID[outliers[id],array],cex=0.8)
text(max(ap)+lambda*0.16,max(b)-mu*0.34,lwd=2,round(xt[xs[outliers[id]],ys[outliers[id]]],2),cex=0.7)
text(max(ap)+lambda*0.16,max(b)-mu*0.36,lwd=2,outliers[id],cex=0.7)


text(max(ap)+lambda*0.08,max(b)-mu*0.1,lwd=2,"x",cex=0.8)
text(max(ap)+lambda*0.08,max(b)-mu*0.14,lwd=2,"y",cex=0.8)
text(max(ap)+lambda*0.08,max(b)-mu*0.18,lwd=2,"R ",cex=0.8)
text(max(ap)+lambda*0.08,max(b)-mu*0.22,lwd=2,"Rb",cex=0.8)
text(max(ap)+lambda*0.08,max(b)-mu*0.26,lwd=2,"Mean",cex=0.8)
text(max(ap)+lambda*0.08,max(b)-mu*0.30,lwd=2,"Code",cex=0.8)
text(max(ap)+lambda*0.08,max(b)-mu*0.34,lwd=2,"Raw",cex=0.7)
text(max(ap)+lambda*0.08,max(b)-mu*0.36,lwd=2,"BeadID",cex=0.7)


if(BLData$G[outliers[id],array]<0){
points(max(ap)+lambda*0.1,max(b)-mu*0.05,cex=2,pch=19,col="yellow")
points(xs[outliers[id]],max(b)-ys[outliers[id]]+min(b),cex=2,pch=21,col="yellow")
}
else{
test<-BLData$G[outliers[id],array]-mean(getProbeIntensities(BLData, BLData$ProbeID[outliers[id],array], array , log=FALSE), na.rm=TRUE)
if(test>0){
points(max(ap)+lambda*0.1,max(b)-mu*0.05,cex=2,pch=19,col="blue")
points(xs[outliers[id]],max(b)-ys[outliers[id]]+min(b),cex=2,pch=21,col="blue")
}
else{
points(max(ap)+lambda*0.1,max(b)-mu*0.05,cex=2,pch=19,col="red")
points(xs[outliers[id]],max(b)-ys[outliers[id]]+min(b),cex=2,pch=21,col="red")
}}


}
else{
text(max(ap)+lambda*0.16,max(b)-mu*0.05,"Spot")
points(max(ap)+lambda*0.1,max(b)-mu*0.05,cex=1.5,pch=3)
points(xs[id],max(b)-ys[id]+min(b),,pch=3,cex=1.5*10/lambda)
points(xs[id],max(b)-ys[id]+min(b),,cex=2,pch=21,col="black")

text(max(ap)+lambda*0.16,max(b)-mu*0.1,lwd=2,xs[id],cex=0.8)
text(max(ap)+lambda*0.16,max(b)-mu*0.14,lwd=2,ys[id],cex=0.8)
text(max(ap)+lambda*0.16,max(b)-mu*0.18,lwd=2,round(BLData$G[id,array],2),cex=0.8)
text(max(ap)+lambda*0.16,max(b)-mu*0.22,lwd=2,round(BLData$Gb[id,array],2),cex=0.8)
text(max(ap)+lambda*0.16,max(b)-mu*0.26,lwd=2,round(mean(getProbeIntensities(BLData,BLData$ProbeID[id,array],array, log=FALSE)),2),cex=0.8)
text(max(ap)+lambda*0.16,max(b)-mu*0.30,lwd=2,BLData$ProbeID[id,array],cex=0.8)
text(max(ap)+lambda*0.16,max(b)-mu*0.34,lwd=2,round(xt[xs[id],ys[id]],2),cex=0.7)
text(max(ap)+lambda*0.16,max(b)-mu*0.36,lwd=2,id,cex=0.7)


text(max(ap)+lambda*0.08,max(b)-mu*0.1,lwd=2,"x",cex=0.7)
text(max(ap)+lambda*0.08,max(b)-mu*0.14,lwd=2,"y",cex=0.7)
text(max(ap)+lambda*0.08,max(b)-mu*0.18,lwd=2,"R",cex=0.7)
text(max(ap)+lambda*0.08,max(b)-mu*0.22,lwd=2,"Rb",cex=0.7)
text(max(ap)+lambda*0.06,max(b)-mu*0.26,lwd=2,"Mean",cex=0.7)
text(max(ap)+lambda*0.08,max(b)-mu*0.30,lwd=2,"code",cex=0.7)
text(max(ap)+lambda*0.08,max(b)-mu*0.34,lwd=2,"TIFFF Int.",cex=0.7)
text(max(ap)+lambda*0.08,max(b)-mu*0.36,lwd=2,"BeadID",cex=0.7)


}

}

"IdentifyIntensity" <-
function(a,ap,b,lambda,mu,loc){
x<-floor(loc$x)
y<-floor(loc$y)
rect(max(ap)+0.01*lambda,min(b)+0.3*mu,max(a)-0.03*lambda,max(b),col="white",border=FALSE)
text(max(ap)+lambda*0.16,max(b)-mu*0.4,lwd=2,round(loc$x,2),cex=0.8)
text(max(ap)+lambda*0.16,max(b)-mu*0.44,lwd=2,round(loc$y,2),cex=0.8)
text(max(ap)+lambda*0.16,max(b)-mu*0.48,lwd=2,round(log2(xt[x,y]),2),cex=0.8)
text(max(ap)+lambda*0.06,max(b)-mu*0.4,lwd=2,"x",cex=0.7)
text(max(ap)+lambda*0.06,max(b)-mu*0.44,lwd=2,"y",cex=0.7)
text(max(ap)+lambda*0.06,max(b)-mu*0.48,lwd=2,"Intensity",cex=0.7)
points(loc$x,loc$y,cex=1,pch=3)

mat<-xt[(x-3):(x+5),(y-6):(y+8)]
for(j in 1:length(mat)){
  if(mat[j] < 1){
 mat[j] <-1
  }}

lambda<-lambda/6
for(i in 1:9){
for(j in 1:15){

colourRamp <- rgb(0,seq(0,1,l=256),0)
rect((max(ap)+(i)*lambda/10+lambda/20),(max(b)-0.05*mu-(j-1)*lambda/10),(max(ap)+(i)*lambda/10+3*lambda/20),(max(b)-(j-1)*lambda/10-0.05*mu+lambda/10),
col=colourRamp[floor(seq(0,256,length=257))[floor(log2(mat[i,15-j])*16)]],border=FALSE)
}
}

}

"FindIlluminaSpotl" <-
function(xt,xs,ys,a,b,...){
indotherspots<-( ( xs>(min(a)-3) ) & ( xs<(max(a)+3) ) & (ys >(min(b)-3) ) & ( ys<(max(b)+3) ) )

    otherspotsX<-(xs[indotherspots])
    otherspotsY<-(ys[indotherspots])
   
    for(j in 1:length(otherspotsX)){
  points(otherspotsX[j],max(b)-otherspotsY[j]+min(b),pch=3,cex=1.5*20/(max(a)-min(a)))  
    }

}




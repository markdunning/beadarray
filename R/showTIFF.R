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
dyellow=(BLData$R[outliers[d],array]<0)
dpos=(BLData$R[outliers[d],array]>0)


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
test[j]=BLData$R[out_blue_red[j],array]-mean(getProbeIntensities(BLData, BLData$ProbeID[out_blue_red[j],array], array, log=FALSE), na.rm=TRUE)
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


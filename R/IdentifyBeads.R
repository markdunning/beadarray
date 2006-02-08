"IdentifyBeads" <-
function(a,ap,b,lambda,mu,id,xt,array,outliers, xs, ys,...){

rect(max(ap)+0.01*lambda,min(b)+0.3*mu,max(a)-0.03*lambda,max(b),col="white",border=FALSE)

if(is.na(match(id,outliers))==FALSE){
outliers[id]<-id
text(max(ap)+lambda*0.16,max(b)-mu*0.05,"Outlier")
text(max(ap)+lambda*0.16,max(b)-mu*0.1,lwd=2,xs[id],cex=0.8)
text(max(ap)+lambda*0.16,max(b)-mu*0.14,lwd=2,ys[id],cex=0.8)
text(max(ap)+lambda*0.16,max(b)-mu*0.18,lwd=2,round(BLData$R[outliers[id],array],2),cex=0.8)
text(max(ap)+lambda*0.16,max(b)-mu*0.22,lwd=2,round(BLData$Rb[outliers[id],array],2),cex=0.8)
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


if(BLData$R[outliers[id],array]<0){
points(max(ap)+lambda*0.1,max(b)-mu*0.05,cex=2,pch=19,col="yellow")
points(xs[outliers[id]],max(b)-ys[outliers[id]]+min(b),cex=2,pch=21,col="yellow")
}
else{
test<-BLData$R[outliers[id],array]-mean(getProbeIntensities(BLData, BLData$ProbeID[outliers[id],array], array , log=FALSE), na.rm=TRUE)
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
text(max(ap)+lambda*0.16,max(b)-mu*0.18,lwd=2,round(BLData$R[id,array],2),cex=0.8)
text(max(ap)+lambda*0.16,max(b)-mu*0.22,lwd=2,round(BLData$Rb[id,array],2),cex=0.8)
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


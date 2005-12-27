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


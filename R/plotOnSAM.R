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


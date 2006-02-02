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


"panel.MA" <-
function(x,y,...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(4,16, -3, 3) )
    points(0.5*(x+y),y-x,pch=16,cex=0.4,col="blue")
    abline(h=c(-1,0,1),lty=c(2,1,2))    
  }


"panel.XY" <-
function(x,y,...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(4,16,4,16) )
    points(x,y,pch=16,cex=0.4)
    abline(0,1,col="orange")    
  }


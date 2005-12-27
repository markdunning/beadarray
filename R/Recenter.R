"Recenter" <-
function(xt,xs,ys,info,outliers,Z){
loc<-locator(1)
x<-loc$x
y<-loc$y
lambda=max(info$a)-min(info$a)
mu=max(info$b)-min(info$b)
anew=(floor(x)-floor(lambda/2)*Z):(floor(x)+floor(lambda/2)*Z)
bnew=(floor(y)-floor(mu/2)*Z):(floor(y)+floor(mu/2)*Z)
showTIFF(xt,xs,ys,anew,bnew, zScale=info$scale, 
                         reverseYaxis = info$Yaxis,
                        width=info$width, height=info$height,
                       resc=info$resc,showSpots=info$showSpots, out=info$out,outliers=outliers, contrast=info$contrast,locate=info$locate,recenter=FALSE)
}


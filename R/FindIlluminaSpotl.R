"FindIlluminaSpotl" <-
function(xt,xs,ys,a,b,...){
indotherspots<-( ( xs>(min(a)-3) ) & ( xs<(max(a)+3) ) & (ys >(min(b)-3) ) & ( ys<(max(b)+3) ) )

    otherspotsX<-(xs[indotherspots])
    otherspotsY<-(ys[indotherspots])
   
    for(j in 1:length(otherspotsX)){
  points(otherspotsX[j],max(b)-otherspotsY[j]+min(b),pch=3,cex=1.5*20/(max(a)-min(a)))  
    }

}


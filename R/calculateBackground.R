"calculateBackground" <-
  function(raw, xs, ys, n=17){
#calculates the background using an n x n square 
    len = length(xs)
    P = vector(length=len)
    n = floor(n / 2)
    for(i in 1:len){
      xc<-xs[i]-floor(xs[i])
      yc<-ys[i]-floor(ys[i])
      dist=vector(length=4)
      dist[1]=xc^2+yc^2
      dist[2]=(xc-1)^2+yc^2
      dist[3]=xc^2+(yc-1)^2
      dist[4]=(xc-1)^2+(yc-1)^2
      cX=c(0,1,0,1)
      cY=c(0,0,1,1)

      newcoord=c((floor(xs[i])+cX[order(dist)[1]]),(floor(ys[i])+cY[order(dist)[1]]))

#Checks to make sure the background never goes out of bounds
      
      if(((newcoord[1]-n) < 0) || ((newcoord[1]+n) > dim(raw)[1]) || ((newcoord[2]-n) < 0) || ((newcoord[2]+n) > dim(raw)[2])){
        P[i] <- mean(P[1]:P[i-1])
        cat("Background value for bead",i,"cannot be calculated, an approximation has been used.\n")
        cat("See manual page for details\n")
      }
      else{
        M<-raw[(newcoord[1]-n):(newcoord[1]+n),(newcoord[2]-n):(newcoord[2]+n)]
        a<-order(M)
        P[i] = round(matrix.mean(M[a[1:5]]),2)
      }
    }
    P
  }

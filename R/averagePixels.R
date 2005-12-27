"averagePixels" <-
  function(xs, ys,raw,k=49777){
    P<-vector(length=k)
    R<-vector(length=k)
    L<-vector(length=k)
    N<-vector(length=k)
    try<-vector(length=k)

    x2 = floor(xs)  
    y2 = floor(ys)

    for(i in 1:k){

      #Checks to make sure the bead isn't extremely close to the edge of the image
      if((x2[i] < 3) || (x2[i] > (dim(raw)[1] - 3)) || (y2[i] < 3) || (y2[i] > (dim(raw)[2] - 3))) {
        L[i] <- NA
        cat("Bead",i,"is too close to the edge of the image to be evaluated and has been ignored.\n")
      }
      else {

        av=vector(length=4)
    
        av[3]=round(matrix.mean(raw[(x2[i]):(x2[i]+2),(y2[i]):(y2[i]+2)]),3)
        av[2]=round(matrix.mean(raw[(x2[i]-1):(x2[i]+1),(y2[i]):(y2[i]+2)]),3)
        av[1]=round(matrix.mean(raw[(x2[i]-1):(x2[i]+1),(y2[i]-1):(y2[i]+1)]),3)
        av[4]=round(matrix.mean(raw[(x2[i]):(x2[i]+2),(y2[i]-1):(y2[i]+1)]),3)

        xc<-xs[i]-floor(xs[i])
        yc<-ys[i]-floor(ys[i])

                                        #w is a set of weights
        w<-c((1-xc)*(1-yc),(1-xc)*yc,xc*yc,xc*(1-yc))

        average=(w[3]*av[3]+w[1]*av[1]+w[4]*av[4]+w[2]*av[2])
        L[i]<-average
      }
    }
    L
  }

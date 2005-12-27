#"sharpen" <-
#	function(xt){

	#Sharpen the intensity of every pixel in the image according to Illumina's algorithm
#
	#Iold<-xt
#	Inew<-Iold

#	a<- length(xt[,1])-1
#	b<- length(xt[1,])-1

#	for(i in 2:a){
#		for(j in 2:b){
#		Inew[i,j] <- Iold[i,j]-0.5*(xt[i+1,j]+xt[i-1,j]+xt[i,j-1]+xt[i,j+1]-4*Iold[i,j]) 
#		}
#	}
#	Inew
#}

sharpen <- function(xt){

  Inew <- xt
  depth <- length(xt[,1])-1

  width <- length(xt[1,])-1
  
  for(i in 2:depth){
    a <- xt[i,2:width]
    b <- xt[i+1,2:width]
    c <- xt[i-1,2:width]
    d <- xt[i,1:(width-1)]
    e <- xt[i,3:(width+1)]

    Inew[i,2:width] <- a - 0.5*(b + c + d + e - (4*a))
  }

  Inew
}

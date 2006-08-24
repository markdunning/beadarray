BLImagePlot <- function(BLData, array = 1, nrow = 18, ncol = 2,
                        low = NULL, high = NULL, ncolors = 123){

  par(mar = c(2,1,1,1), xaxs = "i")
  
#Not needed since the co-ords are automatically scaled to zero now  
#  xs <- floor(BLData$GrnX[,array] - min(BLData$GrnX[,array]))
#  ys <- floor(BLData$GrnY[,array] - min(BLData$GrnY[,array]))

  if (is.character(low)) 
    low <- col2rgb(low)/255
  if (is.character(high)) 
    high <- col2rgb(high)/255
  if (!is.null(low) && is.null(high)) 
    high <- c(1, 1, 1) - low
  if (is.null(low) && !is.null(high)) 
    low <- c(1, 1, 1) - high

  if (is.null(low)) 
    low <- c(1, 1, 1)
  if (is.null(high)) 
    high <- c(0, 1, 0)
  
  col <- rgb(seq(low[1], high[1], len = ncolors), seq(low[2], 
        high[2], len = ncolors), seq(low[3], high[3], len = ncolors))

  xs <- floor(BLData$GrnX[,array])
  ys <- floor(BLData$GrnY[,array])

  xgrid <- floor(seq(0, max(xs), by = max(xs)/ncol))
  ygrid <- floor(seq(0, max(ys), by = max(ys)/nrow))

  imageMatrix <- matrix(ncol = ncol, nrow = nrow)

  for(i in 1:ncol){
    idx = which((xs > xgrid[i]) & (xs < xgrid[i+1]))
    fground = BLData$G[idx,array]
    yvalues = ys[idx]
#    yvalues = BLData$GrnY[idx,array]

    out <- .C("BLImagePlot", length(fground), as.double(fground), as.double(yvalues), as.integer(ygrid),
              result = double(length = nrow), as.integer(nrow), PACKAGE = "beadarray")

    imageMatrix[,i] <- rev(out$result)
  }

  imageMatrix = t((imageMatrix))
  image(x = c(0:ncol), z = imageMatrix,  xaxt = "n", yaxt = "n", col = col)
}

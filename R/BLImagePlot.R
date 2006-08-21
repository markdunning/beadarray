BLImagePlot <- function(BLData, array = 1, rows = 18, columns = 2){
  xs <- floor(BLData$GrnX[,array] - min(BLData$GrnX[,array]))
  ys <- floor(BLData$GrnY[,array] - min(BLData$GrnY[,array]))

  xgrid <- floor(seq(0, max(xs), by = max(xs)/columns))
  ygrid <- floor(seq(0, max(ys), by = max(ys)/rows))

  intensities <- matrix(NA, ncol = columns, nrow = rows)

  for(i in 1:columns){
    for(j in 1:rows){

      indices = which(xs %in% c(xgrid[i]:xgrid[i+1]) & ys %in% c(ygrid[j]:ygrid[j+1]))
      intensities[j,i] = mean(log2(BLData$G[indices,array]),na.rm=TRUE) 
    }
  }
  intensities=t(intensities)
#  print(dim(intensities))
  image(intensities)
}


BLImagePlotC <- function(BLData, array = 1, rows = 18, columns = 2){
  xs <- floor(BLData$GrnX[,array] - min(BLData$GrnX[,array]))
  ys <- floor(BLData$GrnY[,array] - min(BLData$GrnY[,array]))

  xgrid <- floor(seq(0, max(xs), by = max(xs)/columns))
  ygrid <- floor(seq(0, max(ys), by = max(ys)/rows))

  imageMatrix <- matrix(ncol = columns, nrow = rows)

  for(i in 1:columns){
    idx = which((xs > xgrid[i]) & (xs < xgrid[i+1]))
    fground = BLData$G[idx,array]
    yvalues = ys[idx]
#    yvalues = BLData$GrnY[idx,array]

    out <- .C("BLImagePlot", length(fground), as.double(fground), as.double(yvalues), as.integer(ygrid),
              result = double(length = rows), as.integer(rows), PACKAGE = "beadarray")

    imageMatrix[,i] <- rev(out$result)
  }

  imageMatrix = t((imageMatrix))
#  print(dim(imageMatrix))
  image(x = imageMatrix)
}

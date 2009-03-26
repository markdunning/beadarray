HULK <- function(BLData, array, neighbours = NULL, invasions = 20, what = "G") {

  BLD = copyBeadLevelList(BLData)
  arraynms = arrayNames(BLData)

  for(i in array) {
    #print something encouraging
    cat("Array",i,":\n")

    BLD@beadData[[arraynms[i]]]$G = 2^(log2.na(BLD@beadData[[arraynms[i]]]$G) - HULKResids(BLData, i, neighbours, invasions))
  }
  BLD
}

HULKResids <- function(BLData, array, neighbours = NULL, invasions = 20, what = "G") {

  if(is.null(neighbours)) {
      cat("Calculating Neighbourhood\n")
    neighbours <- generateNeighbours(BLData, array)
  }
  residuals <- beadResids(BLData, what = what, array = array)
  weights <- getArrayData(BLData, what = "wts", array = array)
  
  residuals[which((is.na(residuals)) | (weights == 0))] = 0
  
  cat("HULKING\n")
  output <- .C("HULK", as.double(residuals), as.integer(t(neighbours)), as.integer(nrow(neighbours)), as.integer(invasions), results = as.double(residuals))

  output$results
}


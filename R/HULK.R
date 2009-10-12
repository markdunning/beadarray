HULK <- function(BLData, array, neighbours = NULL, invasions = 20) {

  BLD = copyBeadLevelList(BLData)
  arraynms = arrayNames(BLData)

  for(i in array) {
    #print something encouraging
    cat("Array",i,":\n")
    
    if(is.null(neighbours)) {
        cat("Calculating Neighbourhood\n")
        neighbours <- generateNeighbours(BLData, array)
    }
    
    BLD@beadData[[arraynms[i]]]$G = 2^(log2.na(BLD@beadData[[arraynms[i]]]$G) - HULKResids(BLData, i, neighbours, invasions, what = "G"))
    if(BLData@arrayInfo$channels == "two")
        BLD@beadData[[arraynms[i]]]$R = 2^(log2.na(BLD@beadData[[arraynms[i]]]$R) - HULKResids(BLData, i, neighbours, invasions, what = "R"))
  }
  BLD
}

HULKResids <- function(BLData, array, neighbours = NULL, invasions = 20, what = "G") {

  residuals <- beadResids(BLData, what = what, array = array)
  weights <- getArrayData(BLData, what = "wts", array = array)
  
  residuals[which((is.na(residuals)) | (weights == 0))] = 0
  
  if(what == "R")
      cat("HULKING Red Channel\n")
  else
      cat("HULKING Green Channel\n")
      
  output <- .C("HULK", as.double(residuals), as.integer(t(neighbours)), as.integer(nrow(neighbours)), as.integer(invasions), results = as.double(residuals))
  output$results
}



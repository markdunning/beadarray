getProbeIndicesC <- function(BLData, probe, intProbe, index, startSearch = 1){

  if(is.null(BLData$ProbeID)) stop("ProbeID column was not found in BLData object")

  ind <- .C("findIndices", as.integer(probe), intProbe, as.integer(nrow(BLData)), result = integer(length = 25000),
            pos = as.integer(startSearch), PACKAGE="beadarray")

  
  ind2 <- vector()
  i = 1;
  while(ind$result[i] != 0){
    ind2  <- c(ind2, index[ind$result[i]])
    i = i+1;
  }
  list(ind2, ind$pos)
}

getProbeIndicesC <- function(BLData, probe, intProbe){
  ind <- .C("findIndices", as.integer(probe), intProbe, as.integer(nrow(BLData)), result = integer(length = 15000))
  ind2 <- vector()
  i = 1;
  while(ind$result[i] != 0){
    ind2  <- c(ind2, ind$result[i])
    i = i+1;
  }
  ind2
}

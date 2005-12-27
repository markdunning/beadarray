readNonXYData <- function(files, max_length){

##
##files is list of file names to be read
##
##max_length is maximum number of beads on array
##should be 49777 for SAM data
##

  R = matrix(nrow=max_length, ncol=length(files))
  probeID=matrix(nrow=max_length, ncol=length(files), data=0)
  for(i in 1:length(files)){
    r = read.table(files[i], header=T)
    R[1:nrow(r),i] = r[,2]
    probeID[1:nrow(r),i] = r[,1]
  }
  BLData = list()
  BLData$R = R
  BLData$probeID = probeID
  BLData = new("BeadLevelList", BLData)
}

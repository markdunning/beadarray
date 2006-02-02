readProbeInfo <- function(file, beadObject, columns = c("Target", "ProbeID"), header = TRUE, sep = ","){
  info <- read.table(file = file, header = header, sep = sep)
  j <- match(columns, colnames(info),0)
  beadObject$genes <- data.frame(info[,j,drop = FALSE], check.names = FALSE)
  beadObject
}
                     

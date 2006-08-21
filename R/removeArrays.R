removeArrays <- function(BSData, array){
  BSData$G <- BSData$G[,-array]
  BSData$R <- BSData$R[,-array]
  BSData$Gb <- BSData$Gb[,-array]
  BSData$Rb <- BSData$Rb[,-array]
  BSData$BeadStDev <- BSData$BeadStDev[,-array]
  BSData$NoBeads <- BSData$NoBeads[,-array]
  BSData$ProbeID <- BSData$ProbeID[-array]
  BSData$nooutliers <- BSData$nooutliers[,-array]
  if(!is.null(BSData$other)){
    for(i in names(BSData$other)){
      BSData$other[[i]] <- BSData$other[[i]][,-array]
    }
  }
  BSData
}

  

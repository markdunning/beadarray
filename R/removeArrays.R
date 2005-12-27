removeArrays <- function(BSData, array){
  BSData$R <- BSData$R[,-array]
  BSData$G <- BSData$G[,-array]
  BSData$Rb <- BSData$Rb[,-array]
  BSData$Gb <- BSData$Gb[,-array]
  BSData$beadstdev <- BSData$beadstdev[,-array]
  BSData$nobeads <- BSData$nobeads[,-array]
  BSData$probeID <- BSData$probeID[,-array]
  BSData$nooutliers <- BSData$nooutliers[,-array]
  if(!is.null(BSData$other)){
    for(i in names(BSData$other)){
      BSData$other[[i]] <- BSData$other[[i]][,-array]
    }
  }
  BSData
}

  

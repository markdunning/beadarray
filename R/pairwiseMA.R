"pairwiseMA" <-
function(BSData, vec, log= TRUE, labels=c(1:length(vec))){

  if(!class(BSData) == "BeadSummaryList"){
    quit("BeadSummaryList object required!")
  }

  if (log){
    pairs(log2(BSData$R[,vec]),lower.panel=panel.XY,upper.panel=panel.MA, labels=labels)
  }

  else{
    pairs(BSData$R[,vec],lower.panel=panel.XY,upper.panel=panel.MA, labels=labels)
  }
}


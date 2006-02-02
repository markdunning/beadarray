plotMAXY <- function(BSData, vec, log = TRUE, labels=vec,label=FALSE){

  if(!class(BSData) == "BeadSummaryList"){
    quit("BeadSummaryList object required!")
  }

  mat <- matrix(c(0,1,0,0.04, 0,1,0.96,1, 0,0.04,0.04,0.96,
                0.96,1,0.04,0.96, 0.04,0.96,0.04,0.96), byrow = T, ncol= 4)

close.screen(all=TRUE)

  split.screen(mat)

  split.screen(figs = c(length(vec), length(vec)), screen = 5)
  
  for(i in 1:length(vec)){
    for(j in 1:length(vec)){
      screen(((i-1)*length(vec))+j+5)

      par(mar = c(0.3,0.3,0.3,0.3), cex.axis = 0.7)
      if(i == j){
#        plot(0, col.axis = "white", cex = 0, col.lab = "white", tcl = -0, xlab = "", ylab = "")
        plot(0, axes = TRUE, type = "n", tcl = -0, col.axis = "white")
        text(1.0,0, labels = labels[i], cex=2)
      }
      else if(j < i){
        plotXY.beads(BSData, array1 = vec[i], array2 = vec[j], log = log, xaxt = "n", yaxt = "n", label=label)
        if(i == length(vec)){
          axis(1)
          }
        if(j == 1){
          axis(2)
        }
      }
      else{
        plotMA.beads(BSData, array1 = vec[i], array2 = vec[j], log = log, xaxt = "n", yaxt = "n", label=label)
        if(i == 1){
          axis(3)
        }
        if(j == length(vec)){
          axis(4)   
        }
      }
    }
  }
}


"backgroundCorrectBeads" <-
function (BLData, method = "subtract", offset = 0) 
{

if(class(BLData) != "BeadLevelList"){
    stop("BeadLevelList object is required!")
  }

  if(is.null(BLData$Rb)){
    stop("There is no background data for this object")
  }
     
    method <- match.arg(method, c("none", "subtract", "half", 
        "minimum"))
    
    switch(method, subtract = {
        BLData$R <- BLData$R - BLData$Rb
        
    }, half = {
        BLData$R <- pmax(BLData$R - BLData$Rb, 0.5)
        
    }, minimum = {
        BLData$R <- as.matrix(BLData$R - BLData$Rb)
        
for (slide in 1:ncol(BLData$R)) {
            i <- BLData$R[, slide] < 1e-18
            if (any(i, na.rm = TRUE)) {
                m <- min(BLData$R[!i, slide], na.rm = TRUE)
                BLData$R[i, slide] <- m/2
            }
            
        }
    })

    if (offset) {
        BLData$R <- BLData$R + offset
        
    }
    BLData$backgroundCorrected = method

    BLData

}


backgroundCorrect <- 
function (BLData, method = "subtract", offset = 0) 
{

if(class(BLData) != "BeadLevelList"){
    stop("BeadLevelList object is required!")
  }

  if(is.null(BLData@Gb)){
    stop("There is no background data for this object")
  }
     
    method <- match.arg(method, c("none", "subtract", "half", 
        "minimum"))
    
    switch(method, subtract = {
        BLData@G <- BLData@G - BLData@Gb
        
    }, half = {
        BLData@G <- pmax(BLData@G - BLData@Gb, 0.5)
        
    }, minimum = {
        BLData@G <- as.matrix(BLData@G - BLData@Gb)
        
for (slide in 1:ncol(BLData@G)) {
            i <- BLData@G[, slide] < 1e-18
            if (any(i, na.rm = TRUE)) {
                m <- min(BLData@G[!i, slide], na.rm = TRUE)
                BLData@G[i, slide] <- m/2
            }
            
        }
    })

    if (offset) {
        BLData@G <- BLData@G + offset
        
    }
    BLData

}

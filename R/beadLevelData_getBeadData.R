
getBeadData <- function(BLData, what = "Grn", array = 1){
    ##Subset to get all data for the array
    tmp = BLData[[array]]

    m = match.arg(what, colnames(tmp))

    if(is.na(m)) 
        stop("Could not find bead data of type ", what, "\n The available choices are: ", paste(colnames(tmp)), "\n")
        
    return(tmp[,m])
}

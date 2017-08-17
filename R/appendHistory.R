
appendHistory <- function(beadLevelData, text) {
    return(c(beadLevelData@history, paste(Sys.time(), text, sep = ' ')))
}


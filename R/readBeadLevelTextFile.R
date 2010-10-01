numberOfChannels <- function(file, sep = "\t") {
 
    ## Determine the number of channels in a text file based on 
    ## the number of columns
 
    lines <- read.table(file, sep = sep, nrows = 2);
    if(ncol(lines) == 4) ##one channel
        return(1)
    else if(ncol(lines) == 7) ##two channel 
        return(2)
    else ##unexpected number of columns
        return(0)
}

readBeadLevelTextFile <- function(file, sep = "\t") {
    
    ## Read a bead level text file and return a list containing
    ## the contents of the file and how many channels are present
    
    channels <- numberOfChannels(file, sep = sep);
    
    if(channels == 1) 
        data <- matrix(unlist(scan(file, sep = "\t", what = list(integer(), integer(), numeric(), numeric()), skip = 1, quiet = TRUE)), ncol = 4)
    else if (channels == 2)
        data <- matrix(unlist(scan(file, sep = "\t", what = list(integer(), integer(), numeric(), numeric(), integer(), numeric(), numeric()), skip = 1, quiet = TRUE)), ncol = 7)
    else
        stop("Unknown input format!\nExpected 4 columns for single channel data or 7 columns for two channel data\n");
    
    return(data);
}
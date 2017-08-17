numberOfChannels <- function(file, sep = "\t") {

    ## Determine the number of channels in a text file based on 
    ## the number of columns
    
    if(grepl(".bab$", file)) {
        con <- file(file, "rb");
        twoChannel <- (BeadDataPackR:::readHeader( con ))$twoChannel;
        close(con);
        if(twoChannel)
            return(2)
        else
            return(1)
    }
    else {
        lines <- read.table(file, sep = sep, nrows = 2);
        if(ncol(lines) %in% c(4,5)) ##one channel
            return(1)
        else if(ncol(lines) %in% c(7,8)) ##two channel 
            return(2)
        else ##unexpected number of columns
            return(0)
    }
}

numberOfColumns <- function(file, sep = "\t") {
 
    ## Determine the number of channels in a text file based on 
    ## the number of columns
    ## 20-01-11 This should allow the use of weights from Swath data
 
    lines <- read.table(file, sep = sep, nrows = 2);
	if(!all(lines[1,1:4] == c("Code", "Grn", "GrnX", "GrnY"))) 
		return(0)
    ## depending upon whether weights have been included there should be 4/5 cols for single channel
    ## and 7/8 for two channel data
    if(ncol(lines) %in% c(4,5,7,8)) 
        return(ncol(lines))
    else ##unexpected number of columns
        return(0)
}

readBeadLevelTextFile <- function(file, sep = "\t", dec = ".") {
     
    ## Read a bead level text file and return a list containing
    ## the contents of the file and how many channels are present
    ## 20-01-11 Now modified to read a weights column for swath data
     
    columns <- numberOfColumns(file, sep = sep);
     
    if(columns == 4) 
        data <- matrix(unlist(scan(file, sep = sep, what = list(integer(), integer(), numeric(), numeric()), dec = dec, skip = 1, quiet = TRUE)), ncol = 4)
    else if(columns == 5) 
        data <- matrix(unlist(scan(file, sep = sep, what = list(integer(), integer(), numeric(), numeric(), numeric()), dec = dec, skip = 1, quiet = TRUE)), ncol = 5)
    else if(columns == 7) 
        data <- matrix(unlist(scan(file, sep = sep, what = list(integer(), integer(), numeric(), numeric(), integer(), numeric(), numeric()), dec = dec, skip = 1, quiet = TRUE)), ncol = 7)
    else if (columns == 8)
        data <- matrix(unlist(scan(file, sep = sep, what = list(integer(), integer(), numeric(), numeric(), integer(), numeric(), numeric(), numeric()), dec = dec, skip = 1, quiet = TRUE)), ncol = 8)
    else {
        warning("Unknown input format in file", file, "\n  This is probably not a bead-level text file and was ignored\n", call. = FALSE);
		data <- NULL;
	}
    
    return(data);
}

# readBeadLevelTextFile <- function(file, sep = "\t") {
#     
#     ## Read a bead level text file and return a list containing
#     ## the contents of the file and how many channels are present
#     
#     channels <- numberOfChannels(file, sep = sep);
#     
#     if(channels == 1) 
#         data <- matrix(unlist(scan(file, sep = "\t", what = list(integer(), integer(), numeric(), numeric()), skip = 1, quiet = TRUE)), ncol = 4)
#     else if (channels == 2)
#         data <- matrix(unlist(scan(file, sep = "\t", what = list(integer(), integer(), numeric(), numeric(), integer(), numeric(), numeric()), skip = 1, quiet = TRUE)), ncol = 7)
#     else
#         stop("Unknown input format!\nExpected 4 columns for single channel data or 7 columns for two channel data\n");
#     
#     return(data);
# }

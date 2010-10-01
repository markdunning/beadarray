openTIFF <- function(fileName, fileExt)
{
# Provides support for reading the tiff files in compressed formats such as bz2 and gz
# Support for zipped tifs may be added later
    
    con <- switch(fileExt,
                  "tif" = file(fileName, "rb"),
                  "bz2" = bzfile(fileName, "rb"),
                  "gz" = gzfile(fileName, "rb")               
                  )   
    return(con)
}

readTIFF <- function(fileName, path = NULL, verbose = FALSE) {

    # determine the file extension
    fileExt <- unlist(strsplit(fileName, "\\."))
    fileExt <- fileExt[length(fileExt)]
    
    # if there's a PATH argument, use it
    if(!is.null(path)) 
        fileName <- paste(path, fileName, sep = .Platform$file.sep);

    # open the connection
    con <- openTIFF(fileName, fileExt);

    # read the tiff header information
    mode1 <- readBin(con, integer(), 1, size = 1);
    mode2 <- readBin(con, integer(), 1, size = 1);
    version <- readBin(con, integer(), 1, size = 2);
    offset <- readBin(con, integer(), 1, size = 4);

    # we can only use "seek" if the file isn't in a compressed format
    if(fileExt == "tif")
        seek(con, offset)
    else
        tmp <- readBin(con, integer(), size = 1, n = (offset - 8) )

    records <- readBin(con, integer(), 1, size = 2);
    if(verbose) #print information if required
        cat("Mode1:", mode1, "\nMode2:", mode2, "\nVersion:", version, "\nOffset", offset, "\nRecords", records, "\n")
    # loop over each of the records in the header
    for(i in 1:records) {
        tag = readBin(con, integer(), 1, size = 2);
        type = readBin(con, integer(), 1, size = 2);
        length = readBin(con, integer(), 1, size = 4);
        offset <- readBin(con, integer(), 1, size = 4);
    # store the useful information
        if(tag == 256)
            ImageWidth = offset
        if(tag == 257)
            ImageHeight = offset
        if(tag == 273)
            StripOffset = offset
        
        if(verbose)
            cat("Tag:", tag, "Type:", type, "Length:", length, "Offset:", offset, "\n")
    }

    if(fileExt == "tif")
        seek(con, StripOffset)
    else {
        close(con);
        con <- openTIFF(fileName, fileExt);
        tmp <- readBin(con, integer(), size = 1, n = StripOffset)
    }
    
    # read the rest of the image in one go, it's much faster that way!
    data <- readBin(con, integer(), n = ImageWidth * ImageHeight, size = 2, signed = FALSE);
    # rearrange the data into the matrix
    data <- t(matrix(data, ncol = ImageHeight, nrow = ImageWidth));
    close(con);
    return(data);
}

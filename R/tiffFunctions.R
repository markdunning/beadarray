readTIFF <- function(fileName, verbose = FALSE) {
  #open the connection
  con <- file(fileName, "rb");
  #read the tiff header information
  mode1 <- readBin(con, integer(), 1, size = 1);
  mode2 <- readBin(con, integer(), 1, size = 1);
  version <- readBin(con, integer(), 1, size = 2);
  offset <- readBin(con, integer(), 1, size = 4);
  seek(con, offset);
  records <- readBin(con, integer(), 1, size = 2);
  if(verbose) #print information if required
    cat("Mode1:", mode1, "\nMode2:", mode2, "\nVersion:", version, "\nOffset", offset, "\nRecords", records, "\n")
  #loop over each of the records in the header
  for(i in 1:records) {
    tag = readBin(con, integer(), 1, size = 2);
    type = readBin(con, integer(), 1, size = 2);
    length = readBin(con, integer(), 1, size = 4);
    offset <- readBin(con, integer(), 1, size = 4);
    #store the useful information
    if(tag == 256)
      ImageWidth = offset
    if(tag == 257)
      ImageHeight = offset
    if(tag == 273)
      StripOffset = offset

    if(verbose)
      cat("Tag:", tag, "Type:", type, "Length:", length, "Offset:", offset, "\n")
  }
  seek(con, StripOffset);
  #read the rest of the image in one go, it's much faster that way!
  data <- readBin(con, integer(), n = ImageWidth * ImageHeight, size = 2, signed = FALSE);
  #rearrange the data into the matrix
  data <- t(matrix(data, ncol = ImageHeight, nrow = ImageWidth));
  close(con);
  return(data);
}

plotTIFF <- function(tiff, xstart = 0, xend = ncol(tiff)-1, ystart = 0, yend = nrow(tiff)-1, high = "cyan", low = "black", mid = NULL, ncolours = 100, log = TRUE, values = FALSE, textCol = "black", accountForZero = FALSE, ...) {
    
    low = col2rgb(low)/255
    high = col2rgb(high)/255
    if(!is.null(mid))
      mid = col2rgb(mid)/255

    #form the vector of colours across the gradient
    if(is.null(mid))
      colours = rgb(seq(low[1], high[1], len = ncolours), seq(low[2], high[2], len = ncolours), seq(low[3], high[3], len = ncolours))
    else
      colours = c(rgb(seq(low[1], mid[1], len = ncolours/2), seq(low[2], mid[2], len = ncolours/2), seq(low[3], mid[3], len = ncolours/2)), rgb(seq(mid[1], high[1], len = ncolours/2), seq(mid[2], high[2], len = ncolours/2), seq(mid[3], high[3], len = ncolours/2)))
    
        if(log) {
          tiff[which(tiff <= 0)] <- 0.0001
          region <- log2(tiff[(ystart+1):(yend+1), (xstart+1):(xend+1)])
        }
        else 
        region <- tiff[(ystart+1):(yend+1), (xstart+1):(xend+1)]

        if(accountForZero)
                region[which(region <= 0)] <- sort(unique(region))[2]

    colourIndex <- floor((region - min(region)) * ((ncolours-1) / (max(region) - min(region))))+1
    col <- colours[t(colourIndex)]
    
    ##-0.5 adjusts it so the coordinate is the bead centre, not the bottom left corner
    x <- rep(floor(xstart):(floor(xstart)+ncol(region)-1), nrow(region))-0.5
    y <- rep(floor(ystart):floor(ystart+nrow(region)-1), each = ncol(region))-0.5

    plot(0,0, col = "white", xlim = c(min(x), max(x)+1), ylim = c(min(y), max(y)+1), xaxs = "i", yaxs = "i", xlab = "", ylab = "", ...)
    rect(x, y, x+1, y+1, border = NA, col = as.vector(col))
    
    if(values) {
        text(paste(as.vector(t(round(region, 1)))), x = x+0.5, y = y+0.5, col = textCol)
    }  
}
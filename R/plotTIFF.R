plotTIFF <- function(tiff, xrange = c(0, ncol(tiff)-1), yrange = c(0, nrow(tiff)-1), high = "cyan", low = "black", mid = NULL, ncolours = 100, log = TRUE, values = FALSE, textCol = "black", ...) {
    
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
        region <- log2(tiff[(yrange[1]+1):(yrange[2]+1), (xrange[1]+1):(xrange[2]+1)])
    }
    else {
        region <- tiff[(yrange[1]+1):(yrange[2]+1), (xrange[1]+1):(xrange[2]+1)]
    }

    colourIndex <- floor((region - min(region)) * ((ncolours-1) / (max(region) - min(region))))+1
    col <- colours[t(colourIndex)]
    
    ##-0.5 adjusts it so the coordinate is the bead centre, not the bottom left corner
    x <- rep(floor(xrange[1]):(floor(xrange[1])+ncol(region)-1), nrow(region))-0.5
    y <- rep(floor(yrange[1]):floor(yrange[1]+nrow(region)-1), each = ncol(region))-0.5

    plot(0,0, col = "white", xlim = c(min(x), max(x)+1), ylim = c(min(y), max(y)+1), xaxs = "i", yaxs = "i", xlab = "", ylab = "", ...)
    rect(x, y, x+1, y+1, border = NA, col = as.vector(col))
    
    if(values) {
        text(paste(as.vector(t(round(region, 1)))), x = x+0.5, y = y+0.5, col = textCol, ...)
    }  
}

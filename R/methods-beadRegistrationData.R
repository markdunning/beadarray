setMethod("show", signature(object="beadRegistrationData"), function(object) {
    
    cat("Registration Information\n")
    cat("ChipLayout -\n")
    cat(" Two Colour: ", object@layout$twoColor, "\n Number of sections: ", object@layout$nSections,"\n Segments per section: ",  object@layout$nSegs, "\n");
    cat(" Beads per segment: ", nrow(object@coordinateData[[1]]), "\n")
})


setMethod("boxplot",  signature(x="beadRegistrationData"),
    function (x, plotP95 = FALSE, ...) {
        data <- x@registrationData
        if(x@layout$twoColor)
            cols <- rep(rep(c("green", "red"), each = x@layout$nSegs), x@layout$nSections)
        else
            cols <- rep(rep(c("green", "darkgreen"), each = x@layout$nSegs), x@layout$nSections / 2)            

        boxplot(data, col = cols, ylim = c(-5, 10), xaxt = "n", outline = FALSE, ...)
        abline(h = 0, lty = 3)

        if(plotP95) {
#             p95 <- NULL
#             if(x@layout$twoColor) {     
#                     for(i in 1:nrow(x@metrics)) {
#                             p95 <- c(p95, as.integer(x@metrics[i,c("P95Grn", "P95Red")]))
#                     }
#             }
#             else {
#                     p95 <- as.integer(x@metrics[, c("P95Grn")])
#             }
# 
#             p95tmp <- (10 / max(x@p95) ) * p95;
#             lines(x = 1:length(data), y = p95tmp, col = "purple", lwd = 2)           
            p95s <- (10 / max(x@p95) ) * x@p95;
            lines(x = 1:length(data), y = p95s, col = "blue", lwd = 2)
            axis(side = 4, at = c(0, 10), labels = c("0", paste(max(x@p95))))
        }
                   
    }
)

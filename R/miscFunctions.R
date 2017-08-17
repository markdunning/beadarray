log2.na = function (x, ...)
{
    log2(ifelse(x > 0, x, NA), ...)
}

mean.na = function(x) mean(x, na.rm=TRUE)
sd.na = function(x) sd(x, na.rm=TRUE)


####### taken from the sma package #####
plot.smooth.line  <- function(x, M, f = 0.1, ...)
{
#  A <- x
  ind <- !(is.na(x) | is.na(M) | is.infinite(x) | is.infinite(M))
  #lines(lowess(A[ind], M[ind], f = f), ...)
  lines(approx(lowess(x[ind], M[ind], f = f)), ...)  
}

beadarrayUsersGuide <- function(view=TRUE, topic="beadlevel") {
# function modified from limmaUsersGuide() in limma package
        if(!(topic == "beadlevel" || topic == "beadsummary")) {
                cat("\'topic\' must be one of \"beadlevel\" or \"beadsummary\".  Setting topic=\"beadlevel\"\n") 
                topic="beadlevel"
        }
        f = system.file("doc", paste(topic,".pdf",sep=""), package="beadarray")
        if(view) {
                if(.Platform$OS.type == "windows")
                        shell.exec(f)
                else
                        system(paste(Sys.getenv("R_PDFVIEWER"),f,"&"))
        }
        return(f)
}

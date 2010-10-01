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

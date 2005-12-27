#Only use this on matrices that you know aren't going to
#have NA values in.
matrix.mean <- function(x){
  sum(x)/length(x)
}

"read.pgmfile" <-
function (file,...) 
{

fsz <- file.info(file)$size
    con <- file(file, open = "rb")

    pnmhead <- read.pnmhead(con)
 
    ds <- pnmhead$datastart
    seek(con, ds)
    type <- pnmhead$type
    nl <- ifelse(type == "ppm", 3, 1)
    nc <- pnmhead$nc
    nr <- pnmhead$nr
    ncells <- nl * nc * nr
  if (pnmhead$ascii) {
       xx <- scan(con, integer(0), n = ncells)
   }
 #   else {
  #          xx <- readBin(con, "integer", n = ncells, size = 1, 
   #             signed = FALSE)
   # }
    res <- array(xx, dim = c(nc, nr))
   close(con)
res

}


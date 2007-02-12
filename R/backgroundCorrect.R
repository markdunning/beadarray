setGeneric("backgroundCorrect", function(object, method="subtract", offset=0, verbose=FALSE)
   standardGeneric("backgroundCorrect"))

setMethod("backgroundCorrect" ,"BeadLevelList",
          function(object, method = "subtract", offset = 0, verbose=FALSE) {
# Code modified from limma backgroundCorrect function (17 Jan 2007)
    bgc = copyBeadLevelList(object)
    arraynms = arrayNames(object)
    narrays = length(arraynms)
    method = match.arg(method, c("none", "subtract", "half", 
        "minimum", "edwards", "normexp", "rma"))
#    if(is.null(RG$Rb) && is.null(RG$Gb)) method = "none"
    switch(method,
	subtract={
          for(i in 1:narrays) {
            bgc@beadData[[arraynms[i]]]$G = bgc[[arraynms[i]]]$G - bgc[[arraynms[i]]]$Gb
 #           bgc@beadData[[arraynms[i]]]$Gb = NULL
            if(bgc@arrayInfo$channels=="two" | !is.null(bgc[[arraynms[i]]]$Rb)) { # two-colour data
               bgc@beadData[[arraynms[i]]]$R = bgc[[arraynms[i]]]$R - bgc[[arraynms[i]]]$Rb
 #              bgc@beadData[[arraynms[i]]]$Rb = NULL
            }
          }
        },
	half={
          for(i in 1:narrays) {
            bgc@beadData[[arraynms[i]]]$G = pmax(bgc[[arraynms[i]]]$G - bgc[[arraynms[i]]]$Gb, 0.5)
#            bgc@beadData[[arraynms[i]]]$Gb = NULL
            if(bgc@arrayInfo$channels=="two" | !is.null(bgc[[arraynms[i]]]$Rb)) { # two-colour data
               bgc@beadData[[arraynms[i]]]$R = pmax(bgc[[arraynms[i]]]$R - bgc[[arraynms[i]]]$Rb, 0.5)
#               bgc@beadData[[arraynms[i]]]$Rb = NULL
            }
          }
        },
	minimum={
          for(i in 1:narrays) {
            bgc@beadData[[arraynms[i]]]$G = bgc[[arraynms[i]]]$G-bgc[[arraynms[i]]]$Gb
#            bgc@beadData[[arraynms[i]]]$Gb = NULL
            j = bgc[[arraynms[i]]]$G < 1e-18
            if(any(j, na.rm = TRUE)) {
              m = min(bgc[[arraynms[i]]]$G[!j], na.rm = TRUE)
              bgc@beadData[[arraynms[i]]]$G[j] = m/2
            }
            if(bgc@arrayInfo$channels=="two" | !is.null(bgc[[arraynms[i]]]$Rb)) {# two-colour data              bgc@beadData[[arraynms[i]]]$R = bgc[[arraynms[i]]]$R-bgc[[arraynms[i]]]$Rb
#              bgc@beadData[[arraynms[i]]]$Rb = NULL
              j = bgc[[arraynms[i]]]$R < 1e-18
              if(any(j, na.rm = TRUE)) {
                m = min(bgc[[arraynms[i]]]$R[!j], na.rm = TRUE)
                bgc@beadData[[arraynms[i]]]$R[j] = m/2  
              }
            }
          }
        },
	edwards={
#		Log-linear interpolation for dull spots as in Edwards (2003).
#		The threshold values (delta) are chosen such that the number of
#		spots with (0 < R-Rb < delta) is f=10% of the number spots
#		with (R-Rb <= 0) for each channel and array.
#		Note slight change from Edwards (2003).
          for(i in 1:narrays) {
            one = matrix(1,NROW(bgc[[arraynms[i]]]$G),1)
	    delta.vec = function(d, f=0.1) {
	      quantile(d, mean(d<1e-16,na.rm=TRUE)*(1+f), na.rm=TRUE)
	    }
	    sub = as.matrix(bgc[[arraynms[i]]]$G-bgc[[arraynms[i]]]$Gb)
	    delta = one %*% apply(sub, 2, delta.vec)
	    bgc@beadData[[arraynms[i]]]$G = ifelse(sub < delta, delta*exp(1-(bgc[[arraynms[i]]]$Gb+delta)/bgc[[arraynms[i]]]$G), sub)
#            bgc@beadData[[arraynms[i]]]$Gb = NULL
            if(bgc@arrayInfo$channels=="two" | !is.null(bgc[[arraynms[i]]]$Rb)) {# two-colour data
 	      sub = as.matrix(bgc[[arraynms[i]]]$R-bgc[[arraynms[i]]]$Rb)
	      delta = one %*% apply(sub, 2, delta.vec)
	      bgc@beadData[[arraynms[i]]]$R = ifelse(sub < delta, delta*exp(1-(bgc[[arraynms[i]]]$Rb+delta)/bgc[[arraynms[i]]]$R), sub)
#              bgc@beadData[[arraynms[i]]]$Rb = NULL
            }
          }
	},
	normexp={
	  for(i in 1:narrays) {
	    x = bgc[[arraynms[i]]]$G - bgc[[arraynms[i]]]$Gb
	    out = normexp.fit(x)
#	     if(verbose) cat("G: bg.bias=",out$par[1]," bg.sd=",exp(out$par[2])," fg.mean=",exp(out$par[3]),"\n",sep="")
	    bgc@beadData[[arraynms[i]]]$G = normexp.signal(out$par,x)
#            bgc@beadData[[arraynms[i]]]$Gb = NULL
            if(bgc@arrayInfo$channels=="two" | !is.null(bgc[[arraynms[i]]]$Rb)) {# two-colour data
	      x = bgc[[arraynms[i]]]$R-bgc[[arraynms[i]]]$Rb
	      out = normexp.fit(x)
#   	       if(verbose) cat("R: bg.bias=",out$par[1]," bg.log2sd=",out$par[2]," fg.mean=",exp(out$par[3]),"\n",sep="")
	      bgc@beadData[[arraynms[i]]]$R = normexp.signal(out$par,x)
#              bgc@beadData[[arraynms[i]]]$Rb = NULL
            }
	    if(verbose) cat("Corrected array",i,"\n")
	  }
        },
        rma={
	  require("affy")
          for(i in 1:narrays) {      
	    bgc@beadData[[arraynms[i]]]$G = bg.adjust(bgc[[arraynms[i]]]$G-bgc[[arraynms[i]]]$Gb)
#            bgc@beadData[[arraynms[i]]]$Gb = NULL
            if(bgc@arrayInfo$channels=="two" | !is.null(bgc[[arraynms[i]]]$Rb)) {# two-colour data
              bgc@beadData[[arraynms[i]]]$R = bg.adjust(bgc[[arraynms[i]]]$R-bgc[[arraynms[i]]]$Rb)
#              bgc@beadData[[arraynms[i]]]$Rb = NULL
	    }
          }
        },
#        none={
#          for(i in 1:narrays) {
#            bgc@beadData[[arraynms[i]]]$Gb = NULL
#            if(bgc@arrayInfo$channels=="two" | !is.null(bgc[[arraynms[i]]]$Rb)) # two-colour data
#              bgc@beadData[[arraynms[i]]]$Rb = NULL
#          }
#        }
        )    
	if(offset) {
          for(i in 1:narrays) {
            bgc@beadData[[arraynms[i]]]$G = bgc[[arraynms[i]]]$G + offset
            if(bgc@arrayInfo$channels=="two" | !is.null(bgc[[arraynms[i]]]$Rb)) # two-colour data
              bgc@beadData[[arraynms[i]]]$R = bgc[[arraynms[i]]]$R + offset
          }
        }
        bgc
    })

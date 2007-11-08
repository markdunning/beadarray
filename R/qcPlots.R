plotBeadDensities = function(BLData, whatToPlot="G", arrays=NULL, log=TRUE,
   type="l", col=1, xlab="Intensity", ylab="Density", xlim=NULL, ylim=NULL,...) {
   x = y = list()
   arraynms = arrayNames(BLData)
   if(is.null(arrays)) {
     narrays = length(arraynms)
     arrays = 1:narray
   }
   else
     narrays = length(arrays)
   if(length(col)!=narrays*length(whatToPlot))
     col=rep(1, narrays*length(whatToPlot))
   for (i in arrays) {
     for(j in whatToPlot) {
       d = density(getArrayData(BLData, array = i, what = j, log = log), na.rm=TRUE)
       x[[arraynms[i]]][[j]] = d$x
       y[[arraynms[i]]][[j]] = d$y
    }
   }
   if(is.null(xlim))
     xlim = range(x)
   if(is.null(ylim))
     ylim = range(y)
   count = 1
   for(i in arrays) {
     for(j in whatToPlot) {
      if(i==arrays[1] && j==whatToPlot[1])
        plot(x[[arraynms[i]]][[j]], y[[arraynms[i]]][[j]], xlim=xlim, ylim=ylim, col=col[count], xlab=xlab, type=type, ylab=ylab, new=TRUE, ...)
      else
        points(x[[arraynms[i]]][[j]], y[[arraynms[i]]][[j]], col=col[count], type=type, ...)
      count = count+1
    }
 }
}

qcBeadLevel = function(object, whatToPlot="G", RG=FALSE, log=TRUE, nrow= 100, ncol = 100,
                    colDens=1, colBox=0, html=TRUE, fileName="qcsummary.htm", plotdir=NULL, experimentName=NULL, targets=NULL, ...) {

    cat("\nGenerating summary plots\n\n")
    arraynms = arrayNames(object)
    narrays = length(arraynms)
    numplots = length(whatToPlot)
    
    # boxplots of intensities
    if(log)
       ylab = expression(log[2](Intensity))
    else
       ylab = "Intensity"
    for(i in 1:numplots) {
       cat("Boxplots of", whatToPlot[i],"intensities")
       filename = file.path(ifelse(is.null(plotdir), ".", plotdir),
         paste("boxplot", whatToPlot[i], ".png", sep=""))
       png(filename, width=640, height=480)
       boxplotBeads(object, whatToPlot=whatToPlot[i], log=log, outline=FALSE, las=2, main = paste("Boxplots of", whatToPlot[i], "intensities"), ylab=ylab, col=colBox, ...)
       dev.off()
       cat(" .... complete\n")
    }
     
    # density plots of intensities
    if(log)
       xlab = expression(log[2](Intensity))
    else
       xlab = "Intensity"

    for(i in 1:numplots) {
      cat("Density plot of", whatToPlot[i], "intensities")
      filename = file.path(ifelse(is.null(plotdir), ".", plotdir), paste("densities", whatToPlot[i], ".png", sep=""))
      png(filename, width=640, height=480)
      plotBeadDensities(object, whatToPlot=whatToPlot[i], log=log, main=paste("Density plot of", whatToPlot[i], "intensities"), xlab=xlab, col=colDens)
      dev.off()
      cat(" .... complete\n")
      }
      
    # individual image and density plots
    cat("\nGenerating per array plots\n\n")
    for(i in 1:narrays) {
      cat(arraynms[i], ": imageplot ")
      for(j in 1:numplots) {
        filename = file.path(ifelse(is.null(plotdir), ".", plotdir), paste(arraynms[i], ".imageplot", whatToPlot[j], ".png", sep=""))
        png(filename)
        imageplot(object, array=i, whatToPlot=whatToPlot[j], log=log,
                  nrow=nrow, ncol=ncol, ...)
        dev.off()
      }
      cat("density plot ")
      for(j in 1:numplots) {
        filename = file.path(ifelse(is.null(plotdir), ".", plotdir), paste(arraynms[i], ".density", whatToPlot[j], ".png", sep=""))
        png(filename)
        plotBeadDensities(object, array=i, whatToPlot=whatToPlot[j], log=log)
        dev.off()
      }
      if(RG) {
        cat("RG plot ")
        filename = file.path(ifelse(is.null(plotdir), ".", plotdir), paste(arraynms[i], ".RG.png", sep=""))
        png(filename)
        plotRG(object, arrays=i, log=log)
        dev.off()
      }
      cat(" .... complete\n")
    }
      
   if(html) {
    cat("\nCreating an html summary page")
     html = list(NULL)
     
     html[[1]] = "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n<html>\n<head>\n<title>Quality assessment of bead-level data</title>\n<meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\">\n</head>"
     html[[2]] = "<h1>Quality assessment of bead-level data</h1>\n"
     if(!is.null(experimentName))
       html[[2]] = paste(html[[2]],"<h1>",experimentName, "</h1>", sep="")
     html[[3]] = paste("<p>",date(),"</p>\n", sep="")
     html[[4]] = "<h2>Summary plots</h2>\n"
     for(j in 1:numplots) {
       html[[4+j]] = paste("<img src=\"boxplot", whatToPlot[j], ".png\"><br>\n", sep="")
     }
     k = 4+j
     for(j in 1:numplots) {
       html[[k+j]] = paste("<img src=\"densities", whatToPlot[j], ".png\"><br>\n", sep="")
     }
     k = k+j
     html[[k+1]] = "<h2>Individual plots</h2>\n"
     if(is.null(targets)) {
       html[[k+2]] = "<table border=1>\n<tr bgcolor=\"lightgray\"><td><center><strong>Array</strong></center></td><td><center><strong>Plots</strong>,</center/</td></tr>\n"
     }
     else {
       colnms = colnames(targets)
       html[[k+2]] = "<table border=1>\n<tr bgcolor=\"lightgray\">"
       for(i in 1:length(colnms)) {
         html[[k+2]] = paste(html[[k+2]], "<td><center><strong>", colnms[i], "</strong></center></td>", sep="")
        }
        html[[k+2]] = paste(html[[k+2]], "<td><center><strong>Plots</strong></center/</td></tr>\n", sep="")
     }
     for(i in 1:narrays) {
       if(is.null(targets)) {
         html[[k+2+i]] = paste("<tr><td><center>",arraynms[i], "</center></td><td><center>")
       }
       else {
         colnms = colnames(targets)
         html[[k+2+i]] = "<tr><td><center>"
         for(m in 1:length(colnms))
           html[[k+2+i]] = paste(html[[k+2+i]], targets[i,m], "</center></td><td><center>")
       }
       tmp = NULL
       for(j in 1:numplots)
         tmp = paste(tmp, "<a href=\"", arraynms[i], ".imageplot", whatToPlot[j], ".png\">imageplot", whatToPlot[j], "</a>, ", sep="")
       for(j in 1:numplots) {
         if(j==numplots)
           tmp = paste(tmp, "<a href=\"", arraynms[i], ".density", whatToPlot[j], ".png\">density", whatToPlot[j], "</a>", sep="")
         else
           tmp = paste(tmp, "<a href=\"", arraynms[i], ".density", whatToPlot[j], ".png\">density", whatToPlot[j], "</a>, ", sep="")
       }
       if(RG)
         tmp = paste(tmp, ", <a href=\"", arraynms[i], ".RG.png\">RG</a>", sep="")
       tmp = paste(tmp, "</center></td></tr>", sep="")
       html[[k+2+i]] = paste(html[[k+2+i]], tmp, sep="")
       }
     k = k+2+i
     html[[k+1]] = "</table>\n</body>\n</html>"

     filename = file.path(ifelse(is.null(plotdir), ".", plotdir), fileName)
     writeLines(unlist(html), con=filename)    
     cat(" .... complete\n\nOpen", fileName,"to view results.\n")
    }
}

"probeDiagnostics" <-
  function(BLData,probe,array, label=TRUE, log=FALSE, n=3, showOutliers=TRUE,...){

#Makes selection of plots describing the features of a particular probe type and
#their position on the array.

#Check the object is of class BeadLevelList
  if(class(BLData) != "BeadLevelList"){
    stop("BeadLevelList object required!")
  }
    
  #Check that x,y co-ords are available

    if(is.null(BLData$x)){
      stop("X and Y co-ordinates are not present in this data")
    }

    #arrx and arry take all the x and y co-ordinates for every bead in the array
    arrx <- BLData$x[,array]
    arrx<-as.integer(arrx)
    arry <- BLData$y[,array]
    arry<-as.integer(arry)

    par(cex.axis=0.5)
    par(cex.lab=0.3)

    coords = getProbeCoords(BLData, probe, array)

    #xmax and ymax are maximum co-ordinates

    xmax = max(arrx)-min(arrx)
    ymax = max(arry)-min(arry)

    par(mfrow=c(2,2))

#    probePlot(BLData, probe, array, log=log, n=n, label=label, showOutliers=showOutliers,...)
	plotBeadLocations(BLData, probeID=probe, array=array, label=label)
    #plots position of each co-ordinate

    plotCluster(coords)
    #clusters the co-ordinates

    plotDistances(coords)
    #plots distances between co-ordinates

    plotBeadIntensities(BLData, probe, array,label,log=log, n=n,...)
  }



addFeatureData <- function(data,toAdd = c("SYMBOL", "PROBEQUALITY", "CODINGZONE", "PROBESEQUENCE"), annotation = NULL){


  ##If we've supplied a character vector, use these to get data from an annotation package 
  if(is(toAdd, "vector")){


   ##if annotation slot is null, assume it is stored with the object

    if(is.null(annotation)){

      ###should use a getAnnotation function when we have one

      annoName <- annotation(data)

    } else {
      annoName <- annotation
    }

    annoLoaded <- require(paste("illumina", annoName, ".db",sep=""), character.only=TRUE)

    if(annoLoaded){
  
      ##should somehow check that the mapping exists!
    
      mapEnv <-  sapply(paste("illumina", annoName, toAdd,sep=""),as.name)

      IDs <- featureNames(data)

      l <- lapply(mapEnv, function(x) mget(IDs, eval(x), ifnotfound = NA))

    

      newAnno <- data.frame(matrix(unlist(l), nrow = length(IDs), byrow=FALSE))
      rownames(newAnno) <- as.character(IDs)
      colnames(newAnno) <- toAdd

    ###merge the myFeatures data frame

    featureData(data) = new("AnnotatedDataFrame", data=data.frame(merge(fData(data), newAnno, by=0,sort=FALSE), row.names=IDs))

     data

    } else {

      stop("Could not load the annotation package ", paste("illumina", annoName, ".db",sep=""))

    }

  }

  else if (is(toAdd, "data.frame")){

       featureData(data) = new("AnnotatedDataFrame", data=data.frame(merge(fData(data), toAdd, by=0,sort=FALSE), row.names=IDs))

	data

  }

  else stop("The toAdd argument must either be a character vector or a data frame\n")
  

}


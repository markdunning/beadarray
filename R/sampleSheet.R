readSampleSheet <- function(sheet = "sampleSheet.csv"){


  data <- readLines(sheet)

  sampleStart <- grep("\\[Data\\]", data)

  experimentData <- NULL

  if(grep("\\[Header\\]", data[[1]]) != 1 ){


    warning("Expected to see \\[Header\\] in the first line of the sample sheet")
  }
  
   if(length(sampleStart) == 0){

    stop("Expected to see \\[Data\\] somewhere in the sample sheet")
  }

  
   else experimentData[["sampleSheet"]] <- read.csv(sheet, skip = sampleStart)


  for(i in 2:(sampleStart-1)){
    
      field <- strsplit(data[[i]], ",")[[1]][1]
      experimentData[[field]] <- strsplit(data[[i]], ",")[[1]][2]
   
  }



  experimentData
  

}
setGeneric("sampleSheet", function(object) standardGeneric("sampleSheet"))

setGeneric("sampleSheet<-", function(object, value) standardGeneric("sampleSheet<-"))

setMethod("sampleSheet", signature(object = "ExpressionSetIllumina"), function(object) object@experimentData$sampleSheet)
setMethod("sampleSheet", signature(object = "beadLevelData"), function(object) object@experimentData$sampleSheet)

###insertExperimentData option

setReplaceMethod("sampleSheet",
                 signature=signature(
                   object="beadLevelData",
                   value="data.frame"),
                 function(object, value) {
                     object@experimentData$sampleSheet <- value
                     object
                 })

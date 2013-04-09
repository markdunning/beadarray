readSampleSheet <- function(sheet = "sampleSheet.csv"){


  data <- readLines(sheet)

  sampleStart <- grep("\\[Data\\]", data)

  experimentData <- NULL

  if(!grepl("\\[Header\\]", data[[1]])){
      warning("Expected to see \\[Header\\] in the first line of the sample sheet")
  }
  
  if(length(sampleStart) == 0) {
      stop("Expected to see \\[Data\\] somewhere in the sample sheet")
  }
  else {
      experimentData[["sampleSheet"]] <- read.csv(sheet, skip = sampleStart)
  }

  for(i in 2:(sampleStart-1)){
    
      field <- strsplit(data[[i]], ",")[[1]][1]
      experimentData[[field]] <- strsplit(data[[i]], ",")[[1]][2]
   
  }

  experimentData

}
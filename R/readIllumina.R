readIllumina <- function(dir= ".", useImages = FALSE, illuminaAnnotation=NULL, sectionNames = NULL, metricsFile = "Metrics.txt", forceIScan = FALSE, dec = ".", ...) 
{
	
    dir <- normalizePath(dir);
   
    targets <- analyseDirectory(dir = dir, sectionNames = sectionNames, forceIScan = forceIScan, metricsFile = metricsFile)
    metrics = targets$metrics 	
    targets = targets$targets

    ## if there's an .sdf file, read it 
    sdf = NULL
    sdfName = list.files(dir, pattern=".sdf")	
    if(length(sdfName)){ 
        sdf <- simpleXMLparse(readLines(paste(dir, sdfName[1], sep = .Platform$file.sep), warn = FALSE))	
    }
    nSections <- nrow(targets);

    ## report how many channels there are
    nChannels <- numberOfChannels(paste(targets$directory[1], targets$textFile[1], sep = .Platform$file.sep), sep = "\t");
    
    BLData <- new(Class = "beadLevelData");

    BLData = insertSectionData(BLData, what = "Targets", data=targets)
    if(!is.null(metrics)) BLData = insertSectionData(BLData, what="Metrics", data = metrics)


    if(!is.null(sdf)){
        BLData@experimentData$sdfFile <- paste(dir, sdfName, sep= .Platform$file.sep)
        BLData@experimentData$platformClass <- sdf$Class[[1]];
    }



##    BLData@sectionData <- targets[,1:2];        
    nBeads <- vector(length = nSections);

    for(i in 1:nSections) {
        
        message(paste("Processing section ", targets$sectionName[i], sep = ""));
        
        ## check if we've got a .txt or a .bab file here
        if(grepl(".bab", targets$textFile[i])) {
            data <- BeadDataPackR::readCompressedData(inputFile = targets$textFile[i], path = targets$directory[i]);
            ## this will have nondecoded beads so remove them
            data <- data[-which(data[,1] == 0),];
        }
        else {       
            data <- readBeadLevelTextFile(file.path(targets$directory[i], targets$textFile[i]), dec = dec);
			## if the result is NULL this wasn't a read bead-level text file
			## we may need to do some clean up of Metric, sectionData etc
			if(is.null(data)) {
				BLData@sectionData[["Targets"]] <- BLData@sectionData[["Targets"]][-i,];
				if(!is.null(BLData@sectionData[["Metrics"]])) 
					BLData@sectionData[["Metrics"]] <- BLData@sectionData[["Metrics"]][-i,];
				next;
			}
        }
    
        ##record the ProbeIDs, X and Y coords
        BLData <- insertBeadData(BLData, array = i, what = "ProbeID", data = data[,1])
        BLData <- insertBeadData(BLData, array = i, what = "GrnX", data = data[,3])
        BLData <- insertBeadData(BLData, array = i, what = "GrnY", data = data[,4])           

        ## record the number of decoded beads
        nBeads[i] <- nrow(data);
        
        ## read the green images
        if(useImages && !is.null(targets$greenImage[i])) {
            greenImage <- readTIFF(fileName = as.character(targets$greenImage[i]), path = as.character(targets$directory[i]));
            ## there are wrapper functions for these, but using .Call doesn't require
            ## copying the data in the function call
            bg <- .Call("illuminaBackground", greenImage, data[,3:4], PACKAGE = "beadarray")
            greenImage <- .Call("illuminaSharpen", greenImage, PACKAGE = "beadarray");
            fg <- .Call("illuminaForeground", greenImage, data[,3:4], PACKAGE = "beadarray");
            rm(greenImage);
            
            BLData <- insertBeadData(BLData, array = i, what = "Grn", data = fg - bg)
            BLData <- insertBeadData(BLData, array = i, what = "GrnF", data = fg)
            BLData <- insertBeadData(BLData, array = i, what = "GrnB", data = bg)
              
        }
        ## or extract the data from the .txt file
        else {
            BLData <- insertBeadData(BLData, array = i, what = "Grn", data = data[,2])
        }
        
        ## if this is two channel, read the red data too
        if(nChannels == 2) {
        
            BLData <- insertBeadData(BLData, array = i, what = "RedX", data = data[,6])
            BLData <- insertBeadData(BLData, array = i, what = "RedY", data = data[,7])
            
            ## read the images
            if(useImages && !is.null(targets$redImage[i])) {
                image <- readTIFF(fileName = as.character(targets$redImage[i]), path = as.character(targets$directory[i]));
                ## there are wrapper functions for these, but using .Call doesn't require
                ## copying the data in the function call
                bg <- .Call("illuminaBackground", image, data[,6:7], PACKAGE = "beadarray")
                image <- .Call("illuminaSharpen", image, PACKAGE = "beadarray");
                fg <- .Call("illuminaForeground", image, data[,6:7], PACKAGE = "beadarray");
                rm(image);
                
                BLData <- insertBeadData(BLData, array = i, what = "Red", data = fg - bg)
                BLData <- insertBeadData(BLData, array = i, what = "RedF", data = fg)
                BLData <- insertBeadData(BLData, array = i, what = "RedB", data = bg)
                
            }
            ## or extract the data from the .txt file
            else {
                BLData <- insertBeadData(BLData, array = i, what = "Red", data = data[,5])
            }
        }

        ## if the bead-level data contained weights, set them too
        if( ncol(data) %in% c(5,8) )
            BLData <- setWeights(BLData, wts = data[,ncol(data)], array = i);

    }
    
    ## incorporate things into sectionData slot
    ## number of beads
   		
	## sample groupings
	## if we don't have an sdf default to each section being seperate
    sampleGroup <- 1:nrow(targets)

    if(!is.null(sdf)){
		## we need the chip names for the case when multiple chips are being read in one go
		splitSectionNames <- sapply(as.character(BLData@sectionData$Targets$sectionName), strsplit, "_");
		chipNames <- matrix(unlist(splitSectionNames), ncol = length(splitSectionNames[[1]]), byrow = TRUE)[,1]
		## search for the sample group identifiers taken from the sdf
        tmp <- lapply(sdf$SampleLabels$string[[1]], grep, x = targets$sectionName)
        for(i in 1:length(tmp)) {
            for(j in seq(along = tmp[[i]])) {
                sampleGroup[ tmp[[i]][j] ] <- paste(chipNames[ tmp[[i]][j] ], sdf$SampleLabels$string[[1]][i], sep = "_");
            }
        }
    }

    BLData = insertSectionData(BLData, what="SampleGroup", data = data.frame(SampleGroup = sampleGroup))
    BLData = insertSectionData(BLData, what="numBeads", data=data.frame(numBeads = nBeads))

    if(!is.null(illuminaAnnotation)){
            BLData = setAnnotation(BLData, illuminaAnnotation)
    }   	
    else warning("No Illumina annotation was specified. Try setting manually using setAnnotation..\n")

    return(BLData);

}

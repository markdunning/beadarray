analyseDirectory <- function(dir = NULL, twoChannel = NULL, sectionNames = NULL, channel1 = "Grn", channel2 = "Red", txtSuff = ".txt", imgSuff=".tif", locsSuff=".locs", xmlSuff=".xml", metricsFile = "Metrics.txt", ignore = c(".sdf", "fiducial"), forceIScan = FALSE, verbose = FALSE) {

    ## has a directory been specified?
    ## if not, assume working directory
    if(is.null(dir)){
        dir <- getwd()
        if(verbose){cat("No directory supplied, using working directory.\n")}
    }
    
    ## check that directory exists
    if(!file.exists(dir)){ stop("Directory does not exist.\n") }
    
    ## are there any files in the directory?
    fileList<-list.files(dir)
    numFiles<-length(fileList)
    if(numFiles == 0) 
        stop("Directory is empty.\n") 
        
    ## remove any files in the ignore list
    for(i in ignore) {
        if(verbose) { cat("Ignoring files containing \"", i, "\"\n", sep = "") }
        idx <- grep(i, fileList);
        if(length(idx)) 
            fileList <- fileList[-idx];
    }
    
    ## look for a metrics file and read it
    fmet <- grep(metricsFile,fileList)
    if(verbose){ cat("Found",length(fmet),"metrics files in the directory.\n") }

    storemet<-NULL
    for(i in fmet){
        storemet <- rbind(storemet,read.table(file.path(dir,fileList[i],fsep = .Platform$file.sep),header=T,as.is=T,sep="\t"))
    } 
	if(length(fmet)) ## patch from Juerg Straubhaar (BioC-list 2-Aug-2011)
		fileList <- fileList[-fmet];

    ## First lets see if this is BeadScan or iScan data.  
    ## We'll look for the presence of files with "Swath" in them.
    ## we've only seen these in iScan data
    iScan <- any(grepl("Swath", fileList));
    
    ## If this is iScan data, are there two swaths per section?
    ## We'll just check for "Swath2" in file names
    if(iScan) {
        twoSwaths <- any(grepl("Swath2", fileList));
        ## if there aren't two swath files, we can (probably) use the _perBeadFile.txt
        if(!twoSwaths || forceIScan) {
            txtSuff <- "_perBeadFile.txt";
        }
        else {
            swathFileNames <- fileList[grep("Swath", fileList)];
            matchedSwaths <- checkSwathStatus(swathFileNames);
            if(!matchedSwaths) {
                stop("iScan data detected\nThis must first be processed with the function processSwathData()\n");
            }
            ## remove any .txt file that doesn't contain "Swath", i.e beadarray didn't create it
            fileList <- fileList[-which(grepl(".txt", fileList) & !grepl("Swath", fileList))]
        }
    }
    
    ## infer the number of channels 
    fileGreen <- grepl(channel1,fileList)
    fileRed <- grepl(channel2,fileList)
    if(is.null(twoChannel)) {
        if(any(fileRed))
            twoChannel <- TRUE
        else
            twoChannel <- FALSE
    }
    if(verbose){cat("Two Channel:", twoChannel, "\n")}
    
    ## we'll get section names from the .txt file, so see if they exist or are infact .bab
    ## first find the text file   
    #if(!length( grep(txtSuff, fileList) )) {
    #    ## if there aren't any, try and find a bab file instead
    #    if(length( grep(".bab", fileList) )) { 
    #        txtSuff = ".bab";
    #    }
    #    else { ## stop if we still can't find anything
    #        stop(paste("Cannot find a files with extension either", txtSuff, "or .bab"));
    #    }
    #}

	## if there are .bab files, default to those
	## save confusion of targets.txt files etc
	if(length( grep(".bab", fileList) )) {
		txtSuff = ".bab";
	}
	else {
		if(!length( grep(txtSuff, fileList) )) {
			stop(paste("Cannot find any files with extension either", txtSuff, "or .bab"));
		}
	}


    allSections <-  unlist(strsplit(fileList[grep(txtSuff, fileList)], txtSuff))
	
    
    ## if no section names were specified, try to find them all
    if(is.null(sectionNames)) {
        sectionNames <- allSections
    }

   ##Check that if WG-6 data with a _1 or _2 ending are found, the section names we are looking for also have the same extension

    else{
	if(length(grep("_1", allSections)) > 0 | length(grep("_1", allSections)) > 0){
	    if(length(grep("_1", sectionNames)) == 0 & length(grep("_2", sectionNames)) == 0 ){
	      sectionNames <- paste(rep(sectionNames,each=2), c("_1", "_2"),sep="")
	  }
	  
	}
	
	if(length(grep("-Swath1", allSections)) > 0 | length(grep("-Swath2", allSections)) > 0){
	    if(length(grep("-Swath1", sectionNames)) == 0 & length(grep("-Swath2", sectionNames)) == 0 ){
	      sectionNames <- paste(rep(sectionNames,each=2), c("-Swath1", "-Swath2"),sep="")
	  }
	  
	}

    }
      
    
    ## create something to store the results in
    info <- matrix(ncol = 5 + (4^twoChannel), nrow = length(sectionNames));
    colnames(info) <- c("directory", "sectionName", "textFile", "greenImage", "greenLocs", "greenxml", "redImage", "redLocs", "redxml")[1:ncol(info)]
    info[,1] <- rep(dir, nrow(info));
    info[,2] <- sectionNames;
    
    for(i in 1:length(sectionNames)) {
        ## list all the files containing this section name 
        sectionFiles <- fileList[grep(paste(sectionNames[i], "[_.]", sep = ""), fileList)];

        idx <- grep(txtSuff, sectionFiles)
        ## if we cant find a text file for this section, skip it
        if(!length(idx)) {
            message("Data for section ", sectionNames[i], " not found");
            next;
        }
        info[i,3] <- sectionFiles[idx];
            
        ## green images
        idx <- which( grepl(imgSuff, sectionFiles) & grepl(channel1, sectionFiles) )
        if(length(idx)) 
            info[i,4] <- sectionFiles[idx]
        ## green locs    
        idx <- which( grepl(locsSuff, sectionFiles) & grepl(channel1, sectionFiles) )
        if(length(idx)) 
            info[i,5] <- sectionFiles[idx]
        else if (txtSuff == ".bab")
            info[i,5] <- sectionFiles[grep(txtSuff, sectionFiles)]
            
        ## xml file
        if(iScan)
            idx <- grep(xmlSuff, sectionFiles)
        else 
            idx <- which( grepl(xmlSuff, sectionFiles) & grepl(channel1, sectionFiles) )
        if(length(idx))
            info[i,6] <- sectionFiles[idx]
        
        if(twoChannel) {
            idx <- which( grepl(imgSuff, sectionFiles) & grepl(channel2, sectionFiles) )
            if(length(idx))
                info[i,7] <- sectionFiles[idx]
            idx <- which( grepl(locsSuff, sectionFiles) & grepl(channel2, sectionFiles) )            
            if(length(idx))
                info[i,8] <- sectionFiles[idx]
            else if (txtSuff == ".bab")
                info[i,5] <- sectionFiles[grep(txtSuff, sectionFiles)]
                
            ## they seem to only have 1 xml file with iScan, but 2 with BeadScan 
            if(!iScan) {
                idx <- which( grepl(xmlSuff, sectionFiles) & grepl(channel2, sectionFiles) )
                if(length(idx))
                    info[i,9] <- sectionFiles[idx]
            }
            else {
                info <- info[,1:8, drop = FALSE]
            }
                
        }
    }

    ## keep only sections where data was found
    info <- info[which( !is.na(info[,3]) ), ,drop = FALSE]
    if(nrow(info) == 0)
        stop("No data found for the specified sections");
    
    ##match up metrics
    metricsection = "Section"
    metricChip = c("Beadchip", "Matrix");    
    metrow<-rep(NA,length(sectionNames))      
    if(length(fmet)>0){

        if(verbose){cat("Matching up metrics file.\n")}
        if(verbose){cat("Will now look for Section and Chip columns of metrics file.\n")}

        if(!all(colnames(storemet)!=metricsection)){
            if(verbose){cat("Found section column.\n")}

            if(any(colnames(storemet) %in% metricChip)){
                chipIdx <- which(colnames(storemet) %in% metricChip);
                if(verbose){cat("Found chip column.\n")}     

                for(i in 1:nrow(storemet)[1]){

                    matchsect<-grep(paste(storemet[[metricsection]][i], "-?_?", sep = ""), sectionNames)
                    matchchip<-grep(storemet[i,chipIdx],sectionNames)

                    ## this should allow a line in a metrics file to match multiple sections for swath data
                    mymatch<-c(matchchip,matchsect)[duplicated(c(matchchip,matchsect))]
                    for(j in mymatch)
                        metrow[j] <- i   
                }

            }
            
        }
        if(length(metrow)) {
            storemet <- storemet[metrow,,drop = FALSE]
            info <- info[order(metrow), ,drop = FALSE]
        }
        return(list(targets=as.data.frame(info), metrics = as.data.frame(storemet)))
    }
    else {
        return(list(targets = as.data.frame(info)));
    }     
}

expressionQCPipeline = function(BLData, transFun = logGreenChannelTransform, qcDir = "QC", plotType = ".jpeg", horizontal = TRUE,controlProfile=NULL,beadLevelQC = NULL, bashMetrics = NULL, outlierMetrics = NULL, detectionMetrics=NULL,overWrite=FALSE,nSegments=9,outlierFun=illuminaOutlierMethod,tagsToDetect = list(housekeeping = "housekeeping", Biotin = "phage_lambda_genome", Hybridisation = "phage_lambda_genome:high"),zlim=c(5,7),positiveControlTags = c("housekeeping", "phage_lambda_genome"), hybridisationTags =  c("phage_lambda_genome:low", "phage_lambda_genome:med","phage_lambda_genome:high"), negativeTag= "permuted_negative", boxplotFun = logGreenChannelTransform, imageplotFun = logGreenChannelTransform){


an = sectionNames(BLData)

###Use the controlProfile if specified

if(is.null(controlProfile)){

	data(ExpressionControlData)

	controlProfile = makeControlProfile(getAnnotation(BLData))

}

##First step doing per-array plots

dir.create(qcDir, showWarnings=F)
dir.create(paste(qcDir, "/poscont",sep=""), showWarnings=F)
dir.create(paste(qcDir, "/hyb",sep=""), showWarnings=F)
dir.create(paste(qcDir, "/outliers",sep=""), showWarnings=F)
dir.create(paste(qcDir, "/imageplot",sep=""), showWarnings=F)

for(i in 1:length(an)){


	cat("Making per-array plots for section", i, "\n")
	
	##Create a HTML page to write to

	###Plots of control probes

	##Positive controls
	
	cat("Positive controls\n")


	if(plotType == ".jpeg"){

		fname.pc = paste(qcDir, "/poscont/", an[i], ".jpeg",sep="")

		jpeg(fname.pc, width=600, height=400)
	}

	
	else if(plotType == ".pdf"){
		fname.pc = paste(qcDir, "/poscont/", an[i], ".pdf",sep="")	
	 	pdf(fname.pc, width=6, height=4)
	}

	else if(plotType == ".png"){
		fname.pc = paste(qcDir, "/poscont/", an[i], ".png",sep="")
		png(fname.pc, width=600, height=400)
	
	}

	if(file.exists(fname.pc)){
		 if(overWrite){
		poscontPlot(BLData, array=i, controlProfile = controlProfile, positiveControlTags = positiveControlTags, ylim=c(5,16))

		}
		
		else cat("Positive control plot exists. Skipping to next plot\n")	

	}

	else {
	poscontPlot(BLData, array=i, controlProfile = controlProfile, positiveControlTags = positiveControlTags, ylim=c(5,16))

	}
	
	dev.off()
		
	jpeg(paste(qcDir, "/poscont/legend.jpeg", sep=""),width=400, height=400)
	par(mar=c(1,1,1,1))
	plot(1:6, type="n", axes=FALSE, xlab="",ylab="", ylim=c(0,3))

	text(x=c(2,2), y=c(3,2), labels =c("Housekeeping", "Biotin"), adj=0)


	points(x=c(4,4), y=c(3,2), pch = 15, col = c("red", "blue"),cex=2)
	dev.off()
	

	##low/medium/high controls
	
	cat("Hyb controls\n")
	
	if(plotType == ".jpeg"){

		fname.h = paste(qcDir, "/hyb/", an[i], ".jpeg",sep="")

		jpeg(fname.h, width=600, height=400)
	}

	
	else if(plotType == ".pdf"){
		fname.h = paste(qcDir, "/hyb/", an[i], ".pdf",sep="")	
	 	pdf(fname.h, width=6, height=4)
	}

	else if(plotType == ".png"){
		fname.h = paste(qcDir, "/hyb/", an[i], ".png",sep="")
		png(fname.h, width=600, height=400)
	
	}
	

	if(file.exists(fname.h)){
		if(overWrite){
			poscontPlot(BLData, array=i, positiveControlTags =hybridisationTags, controlProfile = controlProfile, colList = NULL, ylim=c(5,16))

		}

		else cat("Hybridisation control plot exists. Skipping to next plot\n")	

	}		

	else{
		poscontPlot(BLData, array=i, positiveControlTags = hybridisationTags, controlProfile = controlProfile, ylim=c(5,16), colList = NULL)

	}
	
	dev.off()

	jpeg(paste(qcDir, "/hyb/legend.jpeg",sep=""), width=400, height=400)
	par(mar=c(1,1,1,1))
	plot(1:6, type="n", axes=FALSE, xlab="",ylab="", ylim=c(0,4))

	text(x=c(2,2,2), y=c(2,3,4), labels =c("Low", "Medium","High"), adj=0)


	points(x=c(4,4,4), y=c(2,3,4), pch = 15, col = c("red", "green","blue"),cex=2)
	dev.off()
	

	cat("Outliers\n")

	if(plotType == ".jpeg") {
		fname.out = paste(qcDir, "/outliers/", an[i], ".jpeg",sep="")
		
		jpeg(fname.out, width=1200, height=300)
	}


	else if(plotType == ".pdf"){
		fname.out = paste(qcDir, "/outliers/", an[i], ".pdf",sep="")
		pdf(fname.out, width=12, height=3)
	}
	else if(plotType == ".png"){
		fname.out = paste(qcDir, "/outliers/", an[i], ".png",sep="")	
 		png(fname.out, width=1200, height=300)
	
	}


	if(file.exists(fname.out)){
		 if(overWrite){
			outlierplot(BLData, array=i, nSegments = nSegments, horizontal = horizontal, outlierFun=outlierFun)
	
		}

		else cat("Outlier plot exists. Skipping to next plot\n")	

	}

	else{
		outlierplot(BLData, array=i, nSegments = nSegments, horizontal = horizontal, outlierFun = outlierFun)
	}

	dev.off()

	cat("imageplot\n")
	
	
	if(plotType == ".jpeg") {
		fname.im = paste(qcDir, "/imageplot/", an[i], ".jpeg",sep="")
		jpeg(fname.im, width=1200, height=300)
	}


	else if(plotType == ".pdf"){
		fname.im = paste(qcDir, "/imageplot/", an[i], ".pdf",sep="")
		pdf(fname.im, width=12, height=3)
	}
	else if(plotType == ".png"){
		fname.im = paste(qcDir, "/imageplot/", an[i], ".png",sep="")	
 		png(fname.im, width=1200, height=300)
	
	}


	if(file.exists(fname.im)){
		if(overWrite){

			imageplot(BLData, array=i, useLocs=TRUE,zlim=zlim, horizontal = horizontal, transFun = imageplotFun)	

		}

		else cat("Positive control plot exists. Skipping to next plot\n")	

	}				
	
	else{
		imageplot(BLData, array=i, useLocs=TRUE,zlim=zlim, horizontal = horizontal, transFun = imageplotFun)	


		
	}

	dev.off()

	###Make a legend for the control plots



	##If array has been BASHed, could do outliers, imageplot and BASHmasked on same plot?


	##Make the HTML page
	outfile = openPage(filename = paste(qcDir, "/",an[i], ".htm", sep=""))
	
	hwrite(paste("Quality assessment for ", an[i]), heading=1,outfile)

	hwrite("Imageplot", heading=2,outfile)

	hwrite("Imageplot created from the log2 transformed green intensiites. White space indicates beads that could not be decoded after array manufacture", outfile)

	hwriteImage(gsub(paste(qcDir, "/", sep=""), "", fname.im),outfile)

	hwrite("Outlier locations", heading=2,outfile)

	hwrite("Locations of beads that are flagged as outliers using Illumina's outlier detection procedure on log2 intensities", outfile)

	hwriteImage(gsub(paste(qcDir, "/", sep=""), "",fname.out),outfile)

	hwrite("Positive Controls", heading=2,outfile)

	hwriteImage(c(gsub(paste(qcDir, "/", sep=""), "",fname.pc), gsub(paste(qcDir, "/", sep=""), "",fname.h)),outfile)
	hwriteImage(c("poscont/legend.jpeg","lmh/legend.jpeg"), outfile)	

	closePage(outfile)


}


outfile = openPage(filename = paste(qcDir, "/Summary.htm", sep=""))


hwrite("Quality assessment summary", heading=1, outfile)


##Create boxplot using defined functions

if(plotType == ".jpeg") jpeg(paste(qcDir, "/Boxplot.jpeg",sep=""), width = 1200, height = 300);hwriteImage("Boxplot.jpeg", outfile)
if(plotType == ".png") png(paste(qcDir, "/Boxplot.jpeg",sep=""), width = 1200, height = 300);hwriteImage("Boxplot.png", outfile)
if(plotType == ".pdf") pdf(paste(qcDir, "/Boxplot.pdf",sep=""), width = 12, height = 3);

boxplot(BLData, transFun = boxplotFun, outline=FALSE)

dev.off()



hwrite("Scan Metrics", heading=2,outfile)

if("Metrics" %in% colnames(BLData@sectionData)) hwrite(BLData@sectionData$Metrics, outfile)

hwrite("Bead-level control summary", heading=2, outfile)

if(is.null(beadLevelQC)){

	cat("Creating probe metrics\n")

	beadLevelQC = makeQCTable(BLData, transFun = transFun, controlProfile = controlProfile)

		
}

else{

	beadLevelQC = getSectionData(BLData, what=beadLevelQC)

}

hwrite(beadLevelQC, outfile)

if(is.null(outlierMetrics)){

	cat("Calculating outlier Metrics\n")

	outlierTable = matrix(nrow = length(an), ncol = nSegments)

	colnames(outlierTable) = paste("Segment", 1:nSegments)
	rownames(outlierTable) = an

	for(i in 1:length(an)){

		outlierTable[i,] = calculateOutlierStats(BLData, transFun = transFun, array=i,nSegments=nSegments, outlierFun=outlierFun)

	}
	
}

else{

	outlierTable = getSectionData(BLData, what=outlierMetrics)

}

hwrite("Outlier Metrics", outfile,heading=2)


hwrite(round(outlierTable,2), outfile)

if(is.null(detectionMetrics)){

	detectionTable = matrix(nrow = length(an), ncol=length(tagsToDetect))
	colnames(detectionTable) = names(tagsToDetect)
	rownames(detectionTable) = an
	for(i in 1:length(an)){
		detectionTable[i,] = controlProbeDetection(BLData, transFun = transFun, array=i, tagsToDetect = tagsToDetect, negativeTag = negativeTag, controlProfile=controlProfile)
	}

}
		
else{
	detectionTable = getSectionData(BLData, what=detectionMetrics)
}

hwrite("Detection Metrics", outfile,heading=2)

hwrite(round(detectionTable, 2),outfile)




closePage(outfile)

if(is.null(bashMetrics)) hwrite("BASH metrics unavailable. Re-run pipeline after BASH\n", outfile)

else{

hwrite(getSectionData(BLData, what=bashMetrics))

}



##Write to csv


write.csv(beadLevelQC, file=paste(qcDir,"/probeMetrics.csv",sep=""), quote=FALSE, row.names=FALSE)

write.csv(BLData@sectionData$Metrics, file=paste(qcDir, "/scanMetrics.csv",sep=""), quote=FALSE, row.names=FALSE)	

write.csv(outlierTable, file=paste(qcDir, "/outlierMetrics.csv",sep=""), quote=FALSE, row.names=FALSE)

write.csv(detectionTable, file=paste(qcDir, "/detectionMetrics.csv",sep=""), quote=FALSE, row.names=FALSE)
	


}

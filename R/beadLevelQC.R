

######Set and get the annotation type of a chip

getAnnotation = function(BLData){

BLData@annotation

}

setAnnotation = function(BLData, aName){

data(ExpressionControlData)

m = match(aName, names(ExpressionControlData))

if(is.na(m)) stop("Annotation must be one of ", paste(names(ExpressionControlData)," "))
else{
#BLData.copy = copyBeadLevelList(BLData)

#BLData.copy@annotation = aName
#BLData.copy
BLData@annotation = aName 
BLData
}
}

getControlAnno = function(BLData){

chipType = getAnnotation(BLData)

if(length(chipType)==0){

stop("No annotation recorded for this BeadLevelList. Try using setAnnotation\n")

}

m = match(chipType, names(ExpressionControlData))

if(is.na(m)) stop("Annotation must be one of ", paste(names(ExpressionControlData)," "))


ExpressionControlData[[match(chipType, names(ExpressionControlData))]]

}

################################################################################
####Calculate detection scores from summarised data


getNegatives = function(BSData){

an = getControlAnno(BSData)

if(is.null(an)) stop("Could not find annotation information for these data")

nIdx = grep("negative", an[,2])

nIDs = match(an[nIdx,1], rownames(exprs(BSData)))

nIDs=nIDs[!is.na(nIDs)]

cat("Number of negative controls found :", length(nIDs),"\n")

if(length(nIDs)!= length(nIdx)) warning("Not all negative controls were found. Please check chip annotation is correct\n")


nData = exprs(BSData)[nIDs, ]

nData


}


calculateDetection = function(BSData){

detScores = matrix(nrow=nrow(exprs(BSData)),ncol=ncol(exprs(BSData)))
 
negControls = getNegatives(BSData)


detect= function(x) 1 - (sum(x>negvals)/(length(negvals)))


for(i in 1:ncol(exprs(BSData))){

negvals = negControls[,i]

M = median(negvals, na.rm=TRUE)
MAD = mad(negvals, na.rm=TRUE)

negvals = negvals[negvals < M+3*MAD]
negvals = negvals[!is.na(negvals)]

detScores[,i] = sapply(exprs(BSData)[,i], detect)

}

detScores
}

###New functionality to calculate QA measures based on bead-level data

###############################################################################
calculateBeadLevelScores=function(BLData,path="QC",log=TRUE,plot=FALSE,replacePlots=TRUE,writeToFile=TRUE){

chipType = getAnnotation(BLData)
data(ExpressionControlData)

if(chipType=="") stop("Could not find annotation type for this chip. Try using setAnnotation\n")


dir.create(path)
dir.create(paste(path, "/outliers",sep=""))
dir.create(paste(path, "/gradients",sep=""))
dir.create(paste(path, "/poscont",sep=""))
#dir.create(paste(QC, "/samplab",sep="")
dir.create(paste(path, "/background",sep=""))
dir.create(paste(path, "/lmh",sep=""))
dir.create(paste(path, "/mismatch",sep=""))



arnames<-arrayNames(BLData)

##need to specify the chipType


outlierMetrics  = matrix(nrow=length(arnames), ncol=36)
controlProbeMetrics = matrix(nrow=length(arnames), ncol=10)

colnames(controlProbeMetrics) = c("HkpDet", "BioDet", "LowDet", "MedDet", "HighDet", "MvsL", "HvsM", "NegDet", "AveNeg", "VarNeg")


colnames(outlierMetrics) = c(paste("Number of beads", 1:9, sep=":"), paste("Number of Illumina outliers", 1:9,sep=":"), paste("Number of low outliers", 1:9,sep=":"), paste("Number of high outliers", 1:9,sep=":"))

rownames(outlierMetrics) = rownames(controlProbeMetrics) = arnames


for(i in 1:length(arnames)){

cat("Producing QC plots for array:", arnames[i], "\n")

outfile = paste(path, "/outliers/",arnames[i],".jpeg",sep="")

cat("Plotting positive controls\n")

outfile = paste(path,"/poscont/",arnames[i],".jpeg",sep="")
controlScores=poscontPlot(BLData,i,plot=plot)

if(plot){


	if(replacePlots){
		jpeg(outfile,width=650,height=375,quality=100)

		controlScores=poscontPlot(BLData,i,plot=plot)

	

		dev.off()
	
	}

	else {
		if(!file.exists(outfile)){
			jpeg(outfile,width=650,height=375,quality=100)
			controlScores=poscontPlot(BLData,i,plot=plot)

			dev.off()
 		}

		else controlScores=poscontPlot(BLData,i,plot=FALSE)
	}
}

else controlScores=poscontPlot(BLData,i,plot=FALSE)





##png(paste("QC/samplab/",arnames[i],".png",sep=""),width=625,height=350,quality=100)
##samplabHV2(BLData,i)
##dev.off()

cat("Plotting hyb controls\n")

outfile = paste(path,"/lmh/",arnames[i],".jpeg",sep="")


if(plot){


	if(replacePlots){
		jpeg(outfile,width=650,height=375,quality=100)

		controlScores=c(controlScores,lmhPlot(BLData,i,plot=plot))


	

		dev.off()
	
	}

	else {
		if(!file.exists(outfile)){
			jpeg(outfile,width=650,height=375,quality=100)
			controlScores=c(controlScores,lmhPlot(BLData,i,plot=plot))


			dev.off()
 		}

		else controlScores=c(controlScores,lmhPlot(BLData,i,plot=FALSE))
	}
}
else controlScores=c(controlScores,lmhPlot(BLData,i,plot=plot))



cat("Plotting mismatch controls\n")

outfile = paste(path, "/mismatch/",arnames[i],".jpeg",sep="")


if(plot){


	if(replacePlots){
		jpeg(outfile,width=650,height=375,quality=100)

		controlScores=c(controlScores,probePairsPlot(BLData,i))


	

		dev.off()
	
	}

	else {
		if(!file.exists(outfile)){
			jpeg(outfile,width=650,height=375,quality=100)
			controlScores=c(controlScores,probePairsPlot(BLData,i))


			dev.off()
 		}

		else controlScores=c(controlScores,probePairsPlot(BLData,i))
	}
}
else controlScores=c(controlScores,probePairsPlot(BLData,i))


cat("Plotting negative controls\n")

outfile = paste(path,"/background/",arnames[i],".jpeg",sep="")

if(plot){


	if(replacePlots){
		jpeg(outfile,width=650,height=375,quality=100)

		controlScores=c(controlScores, backgroundControlPlot(BLData,i,plot=plot))


	

		dev.off()
	
	}

	else {
		if(!file.exists(outfile)){
			jpeg(outfile,width=650,height=375,quality=100)
			controlScores=c(controlScores, backgroundControlPlot(BLData,i,plot=plot))


			dev.off()
 		}

		else controlScores=c(controlScores, backgroundControlPlot(BLData,i,plot=FALSE))
	}
}
else controlScores=c(controlScores, backgroundControlPlot(BLData,i,plot=plot))




controlProbeMetrics[i,] = controlScores


outfile = paste(path,"/outliers/",arnames[i],".jpeg",sep="")

cat("Plotting outlier locations\n")


if(plot){


	if(replacePlots){
		jpeg(outfile,width=1250,height=375,quality=100)

		outlierMetrics[i,]= outlierPlot(BLData,i,log=log,plot=plot)
	

		dev.off()
	
	}

	else {
		if(!file.exists(outfile)){
			jpeg(outfile,width=1250,height=375,quality=100)
			outlierMetrics[i,]= outlierPlot(BLData,i,log=log,plot=plot)
			dev.off()
 		}

		else outlierMetrics[i,]= outlierPlot(BLData,i,log=log,plot=FALSE)

	}
}

else outlierMetrics[i,]= outlierPlot(BLData,i,log=log,plot=FALSE)	



outfile = paste(path,  "/gradients/",arnames[i],".jpeg",sep="")


cat("Plotting chip gradients\n")


if(plot){


	if(replacePlots){
		jpeg(outfile,width=1250,height=375,quality=100)

		gradientPlot(BLData,i,plot=plot)

	

		dev.off()
	
	}

	else {
		if(!file.exists(outfile)){
			jpeg(outfile,width=1250,height=375,quality=100)
			gradientPlot(BLData,i,plot=plot)

			dev.off()
 		}

		else gradientPlot(BLData,i,plot=FALSE)
	}
}

else gradientPlot(BLData,i,plot=FALSE)





if(plot){
cat("Making HTML report\n")

html = list(NULL)
html[[1]] = "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n<html>\n<head>\n<title>Quality assessment of bead-level data</title>\n<meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\">\n</head>"
html[[2]] = paste("<h3>Quality assessment of array ",i,": ", arnames[i],"</h3>", sep = "")
html[[3]] = "<h3>Outliers</h3>\n"
html[[4]] = paste("<img src=\"outliers/",arnames[i],".jpeg\" WIDTH=100%><br>\n",sep="")

html[[5]] = "<h3>Spatial Effects</h3>\n"
html[[6]] = paste("<img src=\"gradients/",arnames[i],".jpeg\" WIDTH=100%><br>\n",sep="")

html[[7]] = "<h3>Control Probes</h3>\n"
html[[8]] = paste("<img src=\"mismatch/",arnames[i],".jpeg\" WIDTH=49%>",sep="")
html[[9]] = paste("<img src=\"background/",arnames[i],".jpeg\" WIDTH=49%><br>\n",sep="")
html[[10]] = paste("<img src=\"lmh/",arnames[i],".jpeg\" WIDTH=49%>",sep="")
html[[11]] = paste("<img src=\"poscont/",arnames[i],".jpeg\" WIDTH=49%><br>\n",sep="")

tnext<-i+1
if(tnext>96){tnext<-tnext-96}
tprev<-i-1
if(tnext<1){tprev<-tprev+96}

html[[12]] = "<center>"
html[[13]] = paste("<a href=\"",arnames[tprev],".htm\">Previous</a> ","<a href=\"",arnames[tnext],".htm\">Next</a>",sep="")
html[[14]] = "<center>"


html[[15]] = "</body>\n</html>"

filename = file.path(path,paste(arnames[i],".htm",sep=""))



writeLines(unlist(html), con=filename)

}

cat("Done QC for array ", arnames[i], "\n")

}

if(writeToFile=="html"){

p=openPage(paste(path,"/qcReport.html",sep=""))

IlluminaOutliers = 100*round(outlierMetrics[,10:18]/outlierMetrics[,1:9],3)
IlluminaOutliers = cbind(IlluminaOutliers, round(100*(apply(outlierMetrics[,10:18],1,sum)/apply(outlierMetrics[,1:9],1,sum)),2))


LowOutliers = 100*round(outlierMetrics[,19:27]/outlierMetrics[,1:9],3)

LowOutliers = cbind(LowOutliers, round(100*(apply(outlierMetrics[,19:27],1,sum)/apply(outlierMetrics[,1:9],1,sum)),2))


HighOutliers = 100*round(outlierMetrics[,28:36]/outlierMetrics[,1:9],3)

HighOutliers = cbind(HighOutliers, round(100*(apply(outlierMetrics[,28:36],1,sum)/apply(outlierMetrics[,1:9],1,sum)),2))


colnames(IlluminaOutliers) = colnames(LowOutliers) = colnames(HighOutliers) = c(paste("Section", 1:9), "Total")


hwrite("Outlier Metrics", p ,heading=1)

hwrite("Outliers found by Illumina", p ,heading=2)

hwrite("Percentage of outliers found using Illumina's default method based on a 3 MAD cut-off from the median of each bead-type.\n 
        Each strip is divided into 9 blocks, and the percentage of beads removed as outliers inside that block is reported",p)

hwrite(IlluminaOutliers, p,col.link=list(paste(path,"/",rownames(IlluminaOutliers),".htm",sep="")))


hwrite("Log-scale outliers below the median", p ,heading=2)

hwrite("Percentages of outliers removed if bead intensities are log-transformed prior to outlier removal. Here percentages are given for beads that are below the median",p)

hwrite(LowOutliers, p,col.link=list(paste(path,"/",rownames(outlierMetrics),".htm",sep="")))

hwrite("Log-scale outliers above the median", p ,heading=2)

hwrite("Percentages of outliers removed if bead intensities are log-transformed prior to outlier removal. Here percentages are given for beads that are above the median",p)

hwrite(HighOutliers, p,col.link=list(paste(path,"/",rownames(outlierMetrics),".htm",sep="")))



hwrite("Control Probe Metrics", p ,heading=2)

hwrite(round(controlProbeMetrics,2), p,col.link=list(paste(path,"/",rownames(controlProbeMetrics),".htm",sep="")))


# scanMetrics = BLData@qcScores$ScanMetrics
scanMetrics = BLData@arrayInfo$qcScores$ScanMetrics

hwrite("Scanner Metrics", p, heading=1)

hwrite(scanMetrics,p,col.link=list(paste(path,"/",rownames(scanMetrics),".htm",sep="")))

closePage(p)

}



#BLData2 = copyBeadLevelList(BLData)
#BLData2@qcScores$OutlierDistribution = outlierMetrics
#BLData2@qcScores$controlProbeScores = controlProbeMetrics
#BLData2
BLData@arrayInfo$qcScores$OutlierDistribution = outlierMetrics
BLData@arrayInfo$qcScores$ControlProbeScores = controlProbeMetrics
BLData



}




probePairsPlot=function(BLData,array=1){

controls = getControlAnno(BLData)


t1<-getArrayData(BLData,what="G",array=array,log=T)
t2<-getArrayData(BLData,what="ProbeID",array=array,log=T)





pms = unique(controls$Array_Address[grep("phage_lambda_genome:pm",controls$Reporter_Group_Identifier)])
mms = unique(controls$Array_Address[grep("phage_lambda_genome:mm2",controls$Reporter_Group_Identifier)])

##which pm probes are found on the array?


pms = pms[pms %in% unique(t2)]
mms = mms[mms %in% unique(t2)]

##Check for unique sequences

uniqueSeqs = unique(controls$Probe_Sequence[match(pms,controls$Array_Address)])

pms = pms[match(uniqueSeqs, controls$Probe_Sequence[match(pms,controls$Array_Address)])]

uniqueSeqs = unique(controls$Probe_Sequence[match(mms,controls$Array_Address)])

mms = mms[match(uniqueSeqs, controls$Probe_Sequence[match(mms,controls$Array_Address)])]


##for each pm probe in turn, find the corresponding mm2

count=1

labs = as.character(pms[1])

pmSeq = controls$Probe_Sequence[match(pms[1], controls$Array_Address)]
mmSeqs = controls$Probe_Sequence[match(mms, controls$Array_Address)]

for(j in 1:length(mms)){

sequenceDiff = sum(unlist(strsplit(as.character(pmSeq),""))!=unlist(strsplit(as.character(mmSeqs[j]),"")))

if(sequenceDiff== 2) mm2 = mms[j]

}


labs = c(labs,as.character(mm2))

newinten = t1[t2 ==  pms[1]]
ys = newinten
xs = rep(count, length(newinten))
cols = rep("red", length(newinten))

count = count +1 

newinten = t1[t2 ==  mm2]
ys = c(ys,newinten)
xs = c(xs,rep(count, length(newinten)))
cols = c(cols,rep("magenta", length(newinten)))

count = count + 1

for(i in 2:length(pms)){

pmSeq = controls$Probe_Sequence[match(pms[i], controls$Array_Address)]
mmSeqs = controls$Probe_Sequence[match(mms, controls$Array_Address)]

for(j in 1:length(mms)){
sequenceDiff = sum(unlist(strsplit(as.character(pmSeq),""))!=unlist(strsplit(as.character(mmSeqs[j]),"")))

if(sequenceDiff== 2) mm2 = mms[j]


}

labs = c(labs, as.character(pms[i]))

newinten = t1[t2 ==  pms[i]]
ys = c(ys, newinten)
xs = c(xs,rep(count, length(newinten)))
cols = c(cols,rep("red", length(newinten)))
count = count +1 

labs = c(labs,as.character(mm2))

newinten = t1[t2 ==  mm2]
ys = c(ys,newinten)
xs = c(xs,rep(count, length(newinten)))
cols = c(cols,rep("magenta", length(newinten)))

count = count + 1

}


plot(0,0,xlim=c(0,count),ylim=c(3,16),type="n",axes=F,ylab="log-intensity",xlab="",main="Low Stringency Hyb Control",xaxs="i")
axis(2)
box()

points(xs+rnorm(length(xs),0,0.02),ys,col=cols,pch=16)
axis(1,labels=labs,at=1:(count-1),las=2)

abline(v=seq(0.5,count,by=2),lty=2)

}



lmhPlot<-function(BLData,array=1,plot=FALSE){

controls = getControlAnno(BLData)


t1<-getArrayData(BLData,what="G",array=array,log=FALSE)
t2<-getArrayData(BLData,what="ProbeID",array=array,log=T)


mrnbg = controls$Array_Address[controls$Reporter_Group_Name=="negative"]

mrnbg = unique(mrnbg[mrnbg %in% unique(t2)])

negvals = exprs(createBeadSummaryData(BLData, array=array, probes=mrnbg))[,1]

M = median(negvals)

MAD = mad(negvals)

outliers = which(negvals>M+3*MAD)

negvals = negvals[-outliers]


detect= function(x) 1 - (sum(x>negvals)/(length(negvals)))

lowID = controls$Array_Address[grep("phage_lambda_genome:low",controls$Reporter_Group_Identifier)]


if(length(is.na(match(lowID, t2)))>0) warning("Could not find all hyb controls on this chip. Please check that it really is a chip of type ", BLData@annotation)

low = t1[t2  %in% lowID]


output= 100*length(which(sapply(low, detect)<0.05))/length(low)

medID = controls$Array_Address[grep("phage_lambda_genome:med",controls$Reporter_Group_Identifier)]

med = t1[t2  %in% medID]

output= c(output,100*length(which(sapply(med, detect)<0.05))/length(med))



highID = controls$Array_Address[grep("phage_lambda_genome:high",controls$Reporter_Group_Identifier)]

high = t1[t2  %in% highID]

output= c(output,100*length(which(sapply(high, detect)<0.05))/length(high))

negvals = low

output= c(output,100*length(which(sapply(med, detect)<0.05))/length(med))

negvals = med

output= c(output,100*length(which(sapply(high, detect)<0.05))/length(high))

names(output) =c("LowDet", "MedDet","HighDet","MvsL","HvsM")


if(plot){
t1=log2(t1)
linepos=0.5

count=1
mrn = controls$Array_Address[grep("phage_lambda_genome:low",controls$Reporter_Group_Identifier)]
mrn = mrn[mrn %in% unique(t2)]

if(length(mrn)==0) stop("Could not find any low hyb controls for this chip. Please check the annotation\n")
newinten = t1[t2 ==  mrn[1]]

labs = mrn[1]

ys = newinten
xs = rep(count, length(newinten))
cols = rep("red", length(newinten))

count = count+1

for(i in 2:length(mrn)){
newinten = t1[t2 == mrn[i]]
ys = c(ys,newinten)
xs =c(xs, rep(count,length(newinten)))
cols = c(cols, rep("red",length(newinten)))
count = count+1
labs = c(labs, mrn[i])
}

##linepos will say where to draw vertical lines

linepos = c(linepos,count-0.5)


mrn = controls$Array_Address[grep("phage_lambda_genome:med",controls$Reporter_Group_Identifier)]
mrn = mrn[mrn %in% unique(t2)]
newinten = t1[t2 ==  mrn[1]]
labs = c(labs, mrn[1])

ys = c(ys, newinten)
xs = c(xs, rep(count, length(newinten)))
cols = c(cols, rep("blue", length(newinten)))

count = count+1

for(i in 2:length(mrn)){
newinten = t1[t2 == mrn[i]]
ys = c(ys,newinten)
xs =c(xs, rep(count,length(newinten)))
cols = c(cols, rep("blue",length(newinten)))
count = count+1
labs = c(labs, mrn[i])

}

linepos=c(linepos,count-0.5)


mrn = controls$Array_Address[grep("phage_lambda_genome:high",controls$Reporter_Group_Identifier)]
mrn = mrn[mrn %in% unique(t2)]
newinten = t1[t2 ==  mrn[1]]
labs = c(labs, mrn[1])

ys = c(ys, newinten)
xs = c(xs, rep(count, length(newinten)))
cols = c(cols, rep("green", length(newinten)))
count = count+1
for(i in 2:length(mrn)){
newinten = t1[t2 == mrn[i]]
ys = c(ys,newinten)
xs =c(xs, rep(count,length(newinten)))
cols = c(cols, rep("green",length(newinten)))
count = count+1
labs = c(labs, mrn[i])

}



plot(0,0,xlim=c(0,count),ylim=c(3,16),type="n",axes=F,ylab="log-intensity",xlab="",main="Hybridisation Controls",xaxs="i")
axis(2)
box()

points(xs+rnorm(length(xs),0,0.02),ys,col=cols,pch=16)
abline(v=linepos,lty=2)
mtext("Low",side=3, at=linepos[1]+0.2)
mtext("Medium",side=3, at=linepos[2]+0.2)
mtext("High",side=3, at=linepos[3]+0.2)
axis(1, labs, at=1:(count-1),las=2)

}

output

}




backgroundControlPlot<-function(BLData,array=1,plot=FALSE){

controls = getControlAnno(BLData)


##Now get 

t1<-getArrayData(BLData,what="G",array=array,log=T)
t2<-getArrayData(BLData,what="ProbeID",array=array,log=T)


mrn = controls$Array_Address[controls$Reporter_Group_Name=="negative"]

mrn=unique(mrn[mrn %in% unique(t2)])
 
###Report the quality scores associated with negatives

negvals = exprs(createBeadSummaryData(BLData, array=array, probes=mrn))[,1]


M = median(negvals)

MAD = mad(negvals)

outliers = which(negvals>M+3*MAD)
negvals = negvals[-outliers]


detect= function(x) 1 - (sum(x>negvals)/(length(negvals)))

totest = t1[t2  %in% mrn]

output= 100*length(which(sapply(negvals, detect)<0.05))/length(negvals)


output= c(output, mean(negvals[-outliers]), var(negvals[-outliers]))

names(output) = c("NegDet","AveNeg","VarNeg")
output


if(plot){

linmat<-matrix(NA,nrow=length(mrn),ncol=3)


t3<-match(t2,mrn)


  for(i in 1:length(mrn)){

  linmat[i,]<-quantile(t1[t2==mrn[i]],c(0.25,0.5,0.75),na.rm=TRUE)
}

linmat<-linmat[order(linmat[,2]),]
  
plot(0,0,type="n",xlim=c(0,7.6),ylim=c(min(linmat,na.rm=T),max(linmat,na.rm=T)),main="negative controls",xlab="Permuted negative controls",ylab="IQR of log-intensity",axes=F)
axis(2)
  box()
  
for(i in 1:length(mrn)){
  lines(c(i/100,i/100),c(linmat[i,1],linmat[i,3]),col="grey50")
}

lines((1:nrow(linmat))/100,linmat[,2])
  
  
}
output
}
  




poscontPlot<-function(BLData,array=1,plot=FALSE){

controls = getControlAnno(BLData)

  t1<-getArrayData(BLData,what="G",array=array,log=FALSE)
  t2<-getArrayData(BLData,what="ProbeID",array=array,log=T)

##this will record where to put vertical lines

linepos=0.5

mrnbg = controls$Array_Address[controls$Reporter_Group_Name=="negative"]

mrnbg=unique(mrnbg[mrnbg %in% unique(t2)])

###Setup average of negatives to be used for detection  

negvals = exprs(createBeadSummaryData(BLData, array=array, probes=mrnbg))[,1]

M = median(negvals)

MAD = mad(negvals)

outliers = which(negvals>M+3*MAD)

negvals = negvals[-outliers]


detect= function(x) 1 - (sum(x>negvals)/(length(negvals)))


bgq<-quantile(t1[match(mrnbg,t2)],c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)

mrnhk = controls$Array_Address[controls$Reporter_Group_Name=="housekeeping"]

mrnhk = mrnhk[mrnhk %in% unique(t2)]

if(length(is.na(match(mrnhk, t2)))>0) warning("Could not find all housekeeping controls on this chip. Please check that it really is a chip of type ", BLData@annotation)


totest = t1[t2  %in% mrnhk]

output= 100*length(which(sapply(totest, detect)<0.05))/length(totest)



mrnb = controls$Array_Address[controls$Reporter_Group_Name=="biotin"]


totest = t1[t2  %in% mrnb]

if(length(is.na(match(mrnb, t2)))>0) warning("Could not find all biotin controls on this chip. Please check that it really is a chip of type ", BLData@annotation)


output= c(output,100*length(which(sapply(totest, detect)<0.05))/length(totest))

names(output) = c("HkpDet","BioDet")

##mrnhs = controls$Array_Address[controls$Reporter_Group_Name=="high_stringency_hyb"]

##totest = t1[t2  %in% mrnhs]

##output= c(output,length(which(sapply(totest, detect)<0.05))/length(totest))

if(plot){

t1=log2(t1)

labs = mrnhk[1]
ys = t1[which(t2 == mrnhk[1])]

if(length(ys)==0) stop("Could not find any Biotin controls on this chip. Please check the annotation\n")

xs = rep(1, length(ys))

cols = rep("red", length(ys))
count =2

if(length(mrnhk)>1){

for(i in 2:length(mrnhk)){

newinten = t1[which(t2 == mrnhk[i])]

ys = c(ys, newinten)
xs =c(xs, rep(count, length(newinten)))
cols = c(cols, rep("red",length(newinten)))
count = count+1
labs = c(labs, mrnhk[i])

}

}


linepos=c(linepos,count-0.5)


mrnb = controls$Array_Address[controls$Reporter_Group_Name=="biotin"]
labs = c(labs, mrnb[1])
mrnb = mrnb[mrnb %in% unique(t2)]



newinten = t1[which(t2 ==  mrnb[1])]
ys = c(ys, newinten)
xs = c(xs, rep(count, length(newinten)))
cols = c(cols, rep("blue", length(newinten)))

count = count+1



for(i in 2:length(mrnb)){
newinten = t1[which(t2 == mrnb[i])]
ys = c(ys, newinten)
xs =c(xs, rep(count, length(newinten)))
cols = c(cols, rep("blue",length(newinten)))
count = count+1

labs = c(labs, mrnb[i])
}

linepos=c(linepos,(count-0.5))

mrnhs = controls$Array_Address[controls$Reporter_Group_Name=="high_stringency_hyb"]
labs = c(labs, mrnhs[1])
mrnhs = mrnhs[mrnhs %in% unique(t2)]
newinten = t1[which(t2 ==  mrnhs[1])]
ys = c(ys, newinten)
xs = c(xs, rep(count, length(newinten)))
cols = c(cols, rep("green", length(newinten)))
count = count +1

if(length(mrnhs)>1){
for(i in 2:length(mrnhs)){
newinten = t1[which(t2 == mrnhs[i])]
ys = c(ys,newinten)
xs =c(xs, rep(count,length(newinten)))
cols = c(cols, rep("green",length(newinten)))
count = count+1
labs = c(labs, mrnhs[i])

}
}

  

plot(0,0,xlim=c(0,count),ylim=c(3,16),type="n",axes=F,ylab="log-intensity",xlab="",main="Positive Controls",xaxs="i")
axis(2)
box()

points(xs+rnorm(length(xs),0,0.02),ys,col=cols,pch=16)
abline(h=bgq,lty=2)
abline(v=linepos,lty=2)
axis(1,labels=labs,at=1:(count-1),las=2)

mtext("Hkeeping",side=3, at=linepos[1]+0.2)
mtext("Bio",side=3, at=linepos[2]+0.2)
mtext("HS",side=3, at=linepos[3]+0.2)

}

output

}



gradientPlot<-function(BLData,array=array,plot=FALSE){

t1<-getArrayData(BLData,what="G",array=array,log=TRUE)
tX<-getArrayData(BLData,what="GrnX",array=array,log=TRUE)
tY<-getArrayData(BLData,what="GrnY",array=array,log=TRUE)

tX = round(100*((tX-min(tX)) / max(tX)))
tY = round(100*((tY-min(tY)) / max(tY)))

###Overall x and y gradients


output=var(sapply(split(t1,tX),mean,na.rm=TRUE),na.rm=TRUE)
output=c(output, var(sapply(split(t1,tY),mean,na.rm=TRUE),na.rm=TRUE))



tYb<-sort(tY)[1+(order((sort(tY))[2:length(tY)]-(sort(tY))[1:(length(tY)-1)],decreasing=T))[1:8]]+sort(tY)[(order((sort(tY))[2:length(tY)]-(sort(tY))[1:(length(tY)-1)],decreasing=T))[1:8]]
tYb<-sort(tYb)/2
gY<-(tY>tYb[1])+(tY>tYb[2])+(tY>tYb[3])+(tY>tYb[4])+(tY>tYb[5])+(tY>tYb[6])+(tY>tYb[7])+(tY>tYb[8])




if(plot){

t1<-getArrayData(BLData,what="residG",array=array,log=TRUE)
tX<-getArrayData(BLData,what="GrnX",array=array,log=TRUE)
tY<-getArrayData(BLData,what="GrnY",array=array,log=TRUE)


par(mar=c(0.1,0.1,3.1,0.1))
layout(matrix(1:27,byrow=T,nrow=3))

#abline(h=0,col="red")
for(i in 0:8){
plot(c(0,100),c(-1,1),type="n",axes=F,main=paste("Y gradient",i+1))
box()
abline(h=0)
plot.smooth.line(100*(tY[gY==i]-min(tY[gY==i],na.rm=T))/(max(tY[gY==i],na.rm=T)-min(tY[gY==i],na.rm=T)),t1[gY==i],col="magenta",lwd=2,f=0.05)
}


for(i in 0:8){
plot(c(0,100),c(-1,1),type="n",axes=F,main=paste("X gradient",i+1))
box()
abline(h=0)
plot.smooth.line(100*(tX[gY==i]-min(tX[gY==i],na.rm=T))/(max(tX[gY==i],na.rm=T)-min(tX[gY==i],na.rm=T)),t1[gY==i],col="blue",lwd=2,f=0.05)
}


for(i in 0:8){
plot(c(0,100),c(-1,1),type="n",axes=F,main=paste("in:out gradient",i+1))
box()
abline(h=0)
mX<-median(tX[gY==i])
mY<-median(tY[gY==i])
tR<-sqrt((tY-mY)^2+(tX-mX)^2)
plot.smooth.line(100*(tR[gY==i]-min(tR[gY==i],na.rm=T))/(max(tR[gY==i],na.rm=T)-min(tR[gY==i],na.rm=T)),t1[gY==i],col="green",lwd=2,f=0.05)
}

}

output

}


outlierPlot<-function(BLData,array=array,log=FALSE,plot=FALSE){
 
tX<-getArrayData(BLData,what="GrnX",array=array,log=TRUE)
tY<-getArrayData(BLData,what="GrnY",array=array,log=TRUE)
resids = beadResids(BLData, array=array,log=log)


tYb<-sort(tY)[1+(order((sort(tY))[2:length(tY)]-(sort(tY))[1:(length(tY)-1)],decreasing=T))[1:8]]+sort(tY)[(order((sort(tY))[2:length(tY)]-(sort(tY))[1:(length(tY)-1)],decreasing=T))[1:8]]
tYb<-sort(tYb)/2
gY<-(tY>tYb[1])+(tY>tYb[2])+(tY>tYb[3])+(tY>tYb[4])+(tY>tYb[5])+(tY>tYb[6])+(tY>tYb[7])+(tY>tYb[8])
  
#fao.ill<-findAllOutliers(BLData,array,log=FALSE,n=3)
#fao.log <-findAllOutliers(BLData,array,log=TRUE,n=3)
fao <-findAllOutliers(BLData,array,log=log,n=3)

#fao.low = intersect(fao.log, which(resids<0))
#fao.high = intersect(fao.log, which(resids>0))
fao.low = intersect(fao, which(resids<0))
fao.high = intersect(fao, which(resids>0))

nbeads = sapply(split(gY,gY),length)
names(nbeads) = paste("Number of beads", 1:9, sep=":")
#outliers.ill=sapply(split(gY[fao.ill],gY[fao.ill]),length)
#names(outliers.ill) = paste("Number of Illumina outliers", 1:9,sep=":")
outliers=sapply(split(gY[fao],gY[fao]),length)
names(outliers) = paste("Number of outliers", 1:9,sep=":")
outliers.low=sapply(split(gY[fao.low],gY[fao.low]),length)
names(outliers.low) = paste("Number of low outliers", 1:9,sep=":")
outliers.high=sapply(split(gY[fao.high],gY[fao.high]),length)
names(outliers.high) = paste("Number of high outliers", 1:9,sep=":")

#output=c(nbeads, outliers.ill,outliers.low,outliers.high)
output=c(nbeads, outliers,outliers.low,outliers.high)

###split the xs into 100 sections
if(plot){


par(mar=c(5.1,0.1,2.1,0.1))
#fao<-findAllOutliers(BLData,array,log=log,n=3)

plot(tY[fao],tX[fao],cex=0.2,pch=16,axes=F,xlab="",ylab="",xaxs="i",yaxs="i")
box()
mtext(paste("no. beads:",length(tX),"no. outliers:",length(fao),"%age outliers:",round(100*length(fao)/length(tX),2)),side=3,line=1)

##fao<-findAllOutliers(BLData,array,log=TRUE,n=4)
##points(tY[fao],tX[fao],cex=0.2,pch=16,col="grey20")
##fao<-findAllOutliers(BLData,array,log=TRUE,n=5)
##points(tY[fao],tX[fao],cex=0.3,pch=16,col="black")

##fao<-findAllOutliers(BLData,array,log=F,n=3)

for(i in 0:8){

  mtext(paste("beads:",length(tX[gY==i])),side=1,line=0,at=median(tY[gY==i]))
  mtext(paste("(",round(length(tX[gY==i])-length(tX)/9),")",sep=""),side=1,line=1,at=median(tY[gY==i]))
  mtext(paste("outliers:",length(fao[gY[fao]==i])),side=1,line=2,at=median(tY[gY==i]))
  mtext(paste("(",round(length(fao[gY[fao]==i])-length(fao)/9),")",sep=""),side=1,line=3,at=median(tY[gY==i]))
  mtext(paste("%age",round(100*length(fao[gY[fao]==i])/length(tX[gY==i]),2)),side=1,line=4,at=median(tY[gY==i]))

  #mtext(paste("(",length(tX),"no. outliers:",length(fao),"%age outliers:",round(100*length(fao)/length(tX),2)),side=1,line=3,at=median(tY[gY==i])
  #mtext(paste("no. beads:",length(tX),"no. outliers:",length(fao),"%age outliers:",round(100*length(fao)/length(tX),2)),side=1,line=4,at=median(tY[gY==i])

  

}

}

output

}



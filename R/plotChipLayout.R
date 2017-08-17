plotChipLayout <-
function(SD,subsC=NULL,sampC=NULL,sectC=NULL,desectC="red",markC="white",main=NULL,samplab=NULL,sectlab=NULL){

## This function takes the sentrix descriptor information collected from the sdf file by simpleXMLparse and plots the chip layout
## Colours and labels are taken from the sdf file unless over-written

## Will only look decent if called to a tall thin plotting device

par(mar=c(4,4,4,4))
plot(c(0,SD$PhysicalProperties$PhysicalWidth[[1]]),c(0,SD$PhysicalProperties$PhysicalHeight[[1]]),type="n",xlab="",ylab="",axes=F,main=main)

# check to see if a user-specified colors are nominated, extract if not
if(is.null(subsC)){subsC<-rgb(SD$UIProperties$SubstrateColor$R[[1]],SD$UIProperties$SubstrateColor$G[[1]],SD$UIProperties$SubstrateColor$B[[1]],maxColorValue=255)}
if(is.null(sampC)){sampC<-rgb(SD$UIProperties$SampleColor$R[[1]],SD$UIProperties$SampleColor$G[[1]],SD$UIProperties$SampleColor$B[[1]],maxColorValue=255)}
if(is.null(sectC)){sectC<-rgb(SD$UIProperties$SectionColor$R[[1]],SD$UIProperties$SectionColor$G[[1]],SD$UIProperties$SectionColor$B[[1]],maxColorValue=255)}

# extract sample dimensions, section dimensions, decode section dimensions and marker dimensions

sampX<-as.numeric(SD$PhysicalProperties$SamplePositions$RelativeXYZPoint$X[[1]])
sampY<-as.numeric(SD$PhysicalProperties$SamplePositions$RelativeXYZPoint$Y[[1]])
sampW<-as.numeric(SD$PhysicalProperties$SampleWidth[[1]])
sampH<-as.numeric(SD$PhysicalProperties$SampleHeight[[1]])
sectX<-as.numeric(SD$PhysicalProperties$SectionPositions$RelativeXYZPoint$X[[1]])
sectY<-as.numeric(SD$PhysicalProperties$SectionPositions$RelativeXYZPoint$Y[[1]])
sectW<-as.numeric(SD$PhysicalProperties$SectionWidth[[1]])
sectH<-as.numeric(SD$PhysicalProperties$SectionHeight[[1]])
desectX<-as.numeric(SD$PhysicalProperties$DecodeSectionPositions$RelativeXYZPoint$X[[1]])
desectY<-as.numeric(SD$PhysicalProperties$DecodeSectionPositions$RelativeXYZPoint$Y[[1]])
desectW<-as.numeric(SD$PhysicalProperties$SectionWidth[[1]])/length(desectX)
desectH<-as.numeric(SD$PhysicalProperties$SectionHeight[[1]])
demarkH<-as.numeric(SD$DecodeMapping$MarkerHeight[[1]])
demarkW<-as.numeric(SD$DecodeMapping$MarkerWidth[[1]])
demarkX<-as.numeric(SD$PhysicalProperties$DecodeMarkerLocations$RelativeXYZPoint$X[[1]])
demarkY<-as.numeric(SD$PhysicalProperties$DecodeMarkerLocations$RelativeXYZPoint$Y[[1]])
anmarkH<-as.numeric(SD$AnalyticalMapping$MarkerHeight[[1]])
anmarkW<-as.numeric(SD$AnalyticalMapping$MarkerWidth[[1]])
anmarkX<-as.numeric(SD$PhysicalProperties$AnalyticalMarkerLocations$RelativeXYZPoint$X[[1]])
anmarkY<-as.numeric(SD$PhysicalProperties$AnalyticalMarkerLocations$RelativeXYZPoint$Y[[1]])


# draw substrate
rect(0,0,SD$PhysicalProperties$PhysicalWidth[[1]],SD$PhysicalProperties$PhysicalHeight[[1]],col=subsC)

# draw samples
for(i in 1:length(sampX)){
	rect(sampX[i]-sampW/2,sampY[i]-sampH/2,sampX[i]+sampW/2,sampY[i]+sampH/2,col=sampC)
	}

#draw sections
for(i in 1:length(sampX)){
	for(j in 1:length(sectX)){
		rect(sampX[i]+sectX[j]-sectW/2,sampY[i]+sectY[j]-sectH/2,sampX[i]+sectX[j]+sectW/2,sampY[i]+sectY[j]+sectH/2,col=sectC)
		}
	}

#draw decode sections
for(i in 1:length(sampX)){
	for(j in 1:length(sectX)){
		for(k in 1:length(desectX)){
			rect(sampX[i]+sectX[j]+desectX[k]-desectW/2,sampY[i]+sectY[j]+desectY[k]-desectH/2,sampX[i]+sectX[j]+desectX[k]+desectW/2,sampY[i]+sectY[j]+desectY[k]+desectH/2,border=desectC)
			}
		}
	}

# draw decode markers
# although described as relative, these seem to be absolute
for(i in 1:length(demarkX)){
	rect(demarkX-demarkW/2,demarkY-demarkH/2,demarkX+demarkW/2,demarkY+demarkH/2,col=markC,border=NA)
	}

# draw analytical markers
# although described as relative, these seem to be absolute
for(i in 1:length(anmarkX)){
	rect(anmarkX-anmarkW/2,anmarkY-anmarkH/2,anmarkX+anmarkW/2,anmarkY+anmarkH/2,col=markC,border=markC)
	}

# check to see if user specified labels supplied, extract if not

if(is.null(samplab)){samplab<-SD$SampleLabels[[1]][[1]]}
if(is.null(sectlab)){sectlab<-SD$SectionLabels[[1]][[1]]}

# draw labels

if((length(unique(sectX))==1)&(length(unique(sampY))==1)){
	axis(side=1,at=sampX,labels=samplab,tick=F)
	mtext("sample",1,0)
	axis(side=2,at=sampY[1]+sectY,labels=sectlab,tick=F,las=1)
	mtext("section",2,2)
	}

## not sure if these work

if((length(unique(sampY))>1)&(length(unique(sampX))==1)){
	axis(side=2,at=sampY,labels=samplab,tick=F,las=1)
	mtext("sample",2,2)
	}

if((length(unique(sampY))>1)&(length(unique(sampX))==2)){
	axis(side=2,at=sampY[sampX==min(sampX)],labels=samplab[sampX==min(sampX)],tick=F,las=1)
	axis(side=4,at=sampY[sampX==max(sampX)],labels=samplab[sampX==max(sampX)],tick=F,las=1)
	}
}


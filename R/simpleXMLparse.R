trimWhiteSpace <- function (x) {
    
    ## code taken from limma
    
    sub("[ \t\n\r]*$", "", sub("^[ \t\n\r]*", "", x))
}


simpleXMLparse <-
function(mysdf){
	
	## This function is designed to take the XML sentrix descriptor file (.sdf) and turn it into a heirarchical list object. I dare say that it may be used for other purposes also
	## It assumes that the sdf file has been read in with a simple readlines command
	## This function requires the limma function trimWhiteSpace	
	
	SentrixDescriptor<-NULL


	## stackdepth maintains a watch on how many layers we are into the XML
	stackdepth<-0

	## Currently if there are more than 12 layers of XML, then the programme cannot cope

	textstring<-rep("",12)
	if(grep("SentrixDescriptor",mysdf[2])==0){stop("FILE NOT IN ANTICIPATED FORMAT")}
	for(line in 3:length(mysdf)){
		#cat(line,"\n")
		linestr<-trimWhiteSpace(mysdf[line])
		while(linestr!=""){

			## The function getnextXMLtoken returns the next token from the sdf file
			## it identifies it as the opening of a new tag, the closing of a tag, or a value to be stored
			## it then returns the remaining string that is still to be processed

			out<-getnextXMLtoken(linestr)
			if(out$stat=="NULL"){
				## if we have a <NAME /> tag, just ignore it and move on
				linestr<-out$rem
				}#end if out$stat = open
			if(out$stat=="OPEN"){
				## if we have encountered a new tag, then we must go a layer deeper into the hierarchy and store the tag name
				stackdepth<-stackdepth+1
				textstring[stackdepth]<-out$val
				linestr<-out$rem
				}#end if out$stat = open
			if(out$stat=="CLOSE"){
				## if we have just closed a tag, then we go a layer back up the hierarchy
				stackdepth<-stackdepth-1
				linestr<-out$rem
				}#end if out$stat = close
			if(out$stat=="VALUE"){	
				## I can't work out how to generalize this, which is why there is a limit of 12 ply depth.
				## This can be easily extended if it is necessary
				## If we have encountered a value to be stored, we must do so in an object with the correct name, alongside anything that is already there

				if(stackdepth==1){SentrixDescriptor[[textstring[1]]]<-list(c((unlist(SentrixDescriptor[[textstring[1]]])),out$val))}

				if(stackdepth==2){SentrixDescriptor[[textstring[1]]][[textstring[2]]]<-list(c((unlist(SentrixDescriptor[[textstring[1]]][[textstring[2]]])),out$val))}

				if(stackdepth==3){SentrixDescriptor[[textstring[1]]][[textstring[2]]][[textstring[3]]]<-list(c((unlist(SentrixDescriptor[[textstring[1]]][[textstring[2]]][[textstring[3]]])),out$val))}

				if(stackdepth==4){SentrixDescriptor[[textstring[1]]][[textstring[2]]][[textstring[3]]][[textstring[4]]]<-list(c((unlist(SentrixDescriptor[[textstring[1]]][[textstring[2]]][[textstring[3]]][[textstring[4]]])),out$val))}

				if(stackdepth==5){SentrixDescriptor[[textstring[1]]][[textstring[2]]][[textstring[3]]][[textstring[4]]][[textstring[5]]]<-list(c((unlist(SentrixDescriptor[[textstring[1]]][[textstring[2]]][[textstring[3]]][[textstring[4]]][[textstring[5]]])),out$val))}

				if(stackdepth==6){SentrixDescriptor[[textstring[1]]][[textstring[2]]][[textstring[3]]][[textstring[4]]][[textstring[5]]][[textstring[6]]]<-list(c((unlist(SentrixDescriptor[[textstring[1]]][[textstring[2]]][[textstring[3]]][[textstring[4]]][[textstring[5]]][[textstring[6]]])),out$val))}

				if(stackdepth==7){SentrixDescriptor[[textstring[1]]][[textstring[2]]][[textstring[3]]][[textstring[4]]][[textstring[5]]][[textstring[6]]][[textstring[7]]]<-list(c((unlist(SentrixDescriptor[[textstring[1]]][[textstring[2]]][[textstring[3]]][[textstring[4]]][[textstring[5]]][[textstring[6]]][[textstring[7]]])),out$val))}

				if(stackdepth==8){SentrixDescriptor[[textstring[1]]][[textstring[2]]][[textstring[3]]][[textstring[4]]][[textstring[5]]][[textstring[6]]][[textstring[7]]][[textstring[8]]]<-list(c((unlist(SentrixDescriptor[[textstring[1]]][[textstring[2]]][[textstring[3]]][[textstring[4]]][[textstring[5]]][[textstring[6]]][[textstring[7]]][[textstring[8]]])),out$val))}
				if(stackdepth==9){SentrixDescriptor[[textstring[1]]][[textstring[2]]][[textstring[3]]][[textstring[4]]][[textstring[5]]][[textstring[6]]][[textstring[7]]][[textstring[8]]][[textstring[9]]]<-list(c((unlist(SentrixDescriptor[[textstring[1]]][[textstring[2]]][[textstring[3]]][[textstring[4]]][[textstring[5]]][[textstring[6]]][[textstring[7]]][[textstring[8]]][[textstring[9]]])),out$val))}

				if(stackdepth==10){SentrixDescriptor[[textstring[1]]][[textstring[2]]][[textstring[3]]][[textstring[4]]][[textstring[5]]][[textstring[6]]][[textstring[7]]][[textstring[8]]][[textstring[9]]][[textstring[10]]]<-list(c((unlist(SentrixDescriptor[[textstring[1]]][[textstring[2]]][[textstring[3]]][[textstring[4]]][[textstring[5]]][[textstring[6]]][[textstring[7]]][[textstring[8]]][[textstring[9]]][[textstring[10]]])),out$val))}

				if(stackdepth==11){SentrixDescriptor[[textstring[1]]][[textstring[2]]][[textstring[3]]][[textstring[4]]][[textstring[5]]][[textstring[6]]][[textstring[7]]][[textstring[8]]][[textstring[9]]][[textstring[10]]][[textstring[11]]]<-list(c((unlist(SentrixDescriptor[[textstring[1]]][[textstring[2]]][[textstring[3]]][[textstring[4]]][[textstring[5]]][[textstring[6]]][[textstring[7]]][[textstring[8]]][[textstring[9]]][[textstring[10]]][[textstring[11]]])),out$val))}

				if(stackdepth==12){SentrixDescriptor[[textstring[1]]][[textstring[2]]][[textstring[3]]][[textstring[4]]][[textstring[5]]][[textstring[6]]][[textstring[7]]][[textstring[8]]][[textstring[9]]][[textstring[10]]][[textstring[11]]][[textstring[12]]]<-list(c((unlist(SentrixDescriptor[[textstring[1]]][[textstring[2]]][[textstring[3]]][[textstring[4]]][[textstring[5]]][[textstring[6]]][[textstring[7]]][[textstring[8]]][[textstring[9]]][[textstring[10]]][[textstring[11]]][[textstring[12]]])),out$val))}
				if(stackdepth>12){stop("stack not deep enough")}
				linestr<-out$rem
				}# end if out$stat = value
				}# while
}
return(SentrixDescriptor)
}


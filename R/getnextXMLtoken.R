getnextXMLtoken <-
function(x){

## This function parses the line of XML passed to it, 
## it extracts and identifies the first XML token in the line
## it returns that token, its class, and the remainder of the line

done<-FALSE
if(nchar(x)>1){
	if(substr(x,1,1)=="<"){
		y<-unlist(strsplit(x,">"))[1]
		if(substr(y,nchar(y),nchar(y))=="/"){
			x<-substr(x,nchar(y)+2,nchar(x))
			type="NULL"
			done<-TRUE
			} # end if substrng ends in a /
		} # end if substring begins with <
	} # end if length of string >0
if(done==FALSE){
	if(nchar(x)>1){
		if(substr(x,1,2)=="</"){
			y<-unlist(strsplit(x,">"))[1]
			x<-substr(x,nchar(y)+2,nchar(x))
			y<-substr(y,3,nchar(y))
			type="CLOSE"
			done<-TRUE
			} # end if substrng begins with a </
		} # end if length of string >0
	} # end if done already
if(done==FALSE){
	if(nchar(x)>1){
		if(substr(x,1,1)=="<"){
			y<-unlist(strsplit(x,">"))[1]
			x<-substr(x,nchar(y)+2,nchar(x))
			y<-substr(y,2,nchar(y))
			type="OPEN"
			done<-TRUE
			} # end if substrng begins with a <
		}# end if length of string >0
	}# end if done already
if(done==FALSE){
	y<-unlist(strsplit(x,"<"))[1]
	x<-substr(x,nchar(y)+1,nchar(x))
	type="VALUE"
	done<-TRUE
	} # end if done already
return(list(rem=x,stat=type,val=y))
}


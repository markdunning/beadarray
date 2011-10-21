

illuminaOutlierMethod= function(inten, probeList,wts=1,n=3)
 {

    probes = sort(unique(probeList[probeList > 0]))

    nasinf = is.na(inten) | !is.finite(inten) | (wts==0)

    inten = inten[!nasinf]
    probeList = probeList[!nasinf]
    nbeads = length(inten)
    start = 0
    foo <- .C("findAllOutliers", as.double(inten), binStatus = integer(length = nbeads), 
        as.integer(probeList), as.integer(probes), as.integer(length(probes)), 
        as.integer(nbeads), as.integer(start), as.double(n), 
        PACKAGE = "beadarray")
    sel = which((probeList > 0) & (foo$binStatus == 0))
    which(!nasinf)[sel]
}

weightsOutlierMethod= function(inten, probeList,wts,n=3)
 {
    probes = sort(unique(probeList[probeList > 0]))
    which(wts==0)
}



noOutlierMethod= function(inten, probeList,wts=1,n=3)
{
integer(0)
}



squeezedVarOutlierMethod<-function (inten, probeList, wts=1, n=3, predictNlim=14){

	# remove na values, infinite values, and non-decoded beads
	
	nasinf = is.na(inten) | !is.finite(inten) | (probeList==0) | (wts==0)
  	inten = inten[!nasinf]
  	probeList = probeList[!nasinf]

	# get initial bead-summary data

  	mysplit<-split(inten,probeList)
  	getmeans<-sapply(mysplit,mean,na.rm=T)
  	getvars<-sapply(mysplit,var,na.rm=T)
  	getN<-sapply(mysplit,length)

	# ensure that we don't extrapolate beyond values based on a decent number of beads

  	usemeans<-getmeans
  	usemeans<-pmax(usemeans,min(getmeans[getN>predictNlim]))
  	usemeans<-pmin(usemeans,max(getmeans[getN>predictNlim]))
  	
  	# model variance in terms of mean
  	
  	myloess<-loess(I(1/getvars[getN>predictNlim])~getmeans[getN>predictNlim])
  	newvars<-predict(myloess,usemeans)
  	
  	# model squared error of variance in terms of mean
  	
  	myloess2<-loess(I((1/getvars[getN>predictNlim]-newvars[getN>predictNlim])^2)~getmeans[getN>predictNlim])
  	newvars2<-predict(myloess2,usemeans)
  	
  	# express priors in terms of precision rather than variance
  	
  	s0<-1/newvars
  	d0<-2*newvars*newvars/newvars2
  	
  	# estimate posterior variance for bead-types
  	
  	shat<-(d0*s0+getN*getvars)/(d0+getN)
  	
  	# construct bead-summary mean and sd for each bead
  	 
  	lsd<-sqrt(shat[match(probeList,unique(probeList))])
  	lmeans<-usemeans[match(probeList,unique(probeList))]
  	
  	# generate list of outliers
  	
  	sel = which(abs((inten-lmeans)/lsd)>n)
  	which(!nasinf)[sel]
}      

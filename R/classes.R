setClass("LargeDataObject")

setClass("BeadSummaryList",
representation("list")
)

setClass("BeadLevelList",
         representation("list")
         )

setIs("BeadLevelList", "LargeDataObject")

setIs("BeadSummaryList", "LargeDataObject")

dim.BeadSummaryList <- function(x)  if(is.null(x[[1]])) c(0,0) else dim(as.matrix(x[[1]]))
length.BeadSummaryList <- length.BeadLevelList <- function(x) prod(dim(x))
dimnames.BeadSummaryList <- function(x) dimnames(x[[1]])

dim.BeadLevelList <- function(x) if(is.null(x$R)) c(0,0) else dim(as.matrix(x$R))
dimnames.BeadLevelList <- function(x) dimnames(x$R)


#allows the subsetting of the BeadSummaryList object.  
assign("[.BeadSummaryList",
function(object, i, j, ...) {

	if (nargs() != 3) stop("Two subscripts required",call.=FALSE)
	oc <- names(object$other)        
	if (missing(i))
		if (missing(j))
			return(object)
		else {
                  object$R <- object$R[,j, drop=FALSE]
                  object$G <- object$G[,j, drop=FALSE]
                  object$Rb <- object$Rb[,j, drop=FALSE]
                  object$Gb <- object$Gb[,j, drop=FALSE]
                  object$BeadStDev <- object$BeadStDev[,j, drop=FALSE]
                  object$Nobeads <- object$Nobeads[,j, drop=FALSE]
			object$Detection <- object$Detection[,j, drop=FALSE]
                  object$ProbeID <- object$ProbeID
                  object$nooutliers <- object$nooutliers[,j, drop=FALSE]
                  if(!is.null(oc)) for(k in oc) object$other[[k]] <- object$other[[k]][,j,drop=FALSE]
                }
	else
		if (missing(j)) {
                  object$R <- object$R[i,, drop=FALSE]
                  object$G <- object$G[i,, drop=FALSE]
                  object$Rb <- object$Rb[i,, drop=FALSE]
                  object$Gb <- object$Gb[i,, drop=FALSE]
                  object$BeadStDev <- object$BeadStDev[i,, drop=FALSE]
                  object$Nobeads <- object$Nobeads[i,, drop=FALSE]
                  object$ProbeID <- object$ProbeID[i]
                  object$nooutliers <- object$nooutliers[i,, drop=FALSE]
                  if(!is.null(oc)) for(k in oc) object$other[[k]] <- object$other[[k]][i,,drop=FALSE]
		} else {
                  object$R <- object$R[i,j, drop=FALSE]
                  object$G <- object$G[i,j, drop=FALSE]
                  object$Rb <- object$Rb[i,j, drop=FALSE]
                  object$Gb <- object$Gb[i,j, drop=FALSE]
                  object$BeadStDev <- object$BeadStDev[i,j, drop=FALSE]
                  object$Nobeads <- object$Nobeads[i,j, drop=FALSE]
                  object$ProbeID <- object$ProbeID[i]
                  object$nooutliers <- object$nooutliers[i,j, drop=FALSE]
                  if(!is.null(oc)) for(k in oc) object$other[[k]] <- object$other[[k]][i,j,drop=FALSE]
		}
	object
})

#allows the subsetting of the BeadLevelList object.
#Modified from limma source code.
assign("[.BeadLevelList",
function(object, i, j, ...) {
	oc <- names(object$other)   
	if (nargs() != 3) stop("Two subscripts required",call.=FALSE)
	if (missing(i))		if (missing(j))
			return(object)
		else {
                  object$R <- object$R[,j, drop=FALSE]
                  object$G <- object$G[,j, drop=FALSE]
                  object$Rb <- object$Rb[,j, drop=FALSE]
                  object$Gb <- object$Gb[,j, drop=FALSE]
                  object$x <- object$x[,j, drop=FALSE]
                  object$y <- object$y[,j, drop=FALSE]
                  object$ProbeID <- object$ProbeID[,j, drop=FALSE]
                  object$targets <- object$targets[j,,drop=FALSE]
                  object$SAMPLE <- object$SAMPLE[j,,drop=FALSE]       
                  if(!is.null(oc)) for(k in oc) object$other[[k]] <- object$other[[k]][,j,drop=FALSE]
		}
	else
		if (missing(j)) {
                  object$R <- object$R[i,, drop=FALSE]
                  object$G <- object$G[i,, drop=FALSE]
                  object$Rb <- object$Rb[i,, drop=FALSE]
                  object$Gb <- object$Gb[i,, drop=FALSE]
                  object$x <- object$x[i,, drop=FALSE]
                  object$y <- object$y[i,, drop=FALSE]
                  object$ProbeID <- object$ProbeID[i,, drop=FALSE]
                  object$SAMPLE <- object$SAMPLE[i,,drop=FALSE]
                  if(!is.null(oc)) for(k in oc) object$other[[k]] <- object$other[[k]][i,,drop=FALSE]
   		} else {
                  object$R <- object$R[i,j, drop=FALSE]
                  object$G <- object$G[i,j, drop=FALSE]
                  object$Rb <- object$Rb[i,j, drop=FALSE]
                  object$Gb <- object$Gb[i,j, drop=FALSE]
                  object$x <- object$x[i,j, drop=FALSE]
                  object$y <- object$y[i,j, drop=FALSE]
                  object$ProbeID <- object$ProbeID[i,j, drop=FALSE]
                  object$targets <- object$targets[j,,drop=FALSE]
                  object$SAMPLE <- object$SAMPLE[j,,drop=FALSE]  
                  if(!is.null(oc)) for(k in oc) object$other[[k]] <- object$other[[k]][i,j,drop=FALSE]
                }
	object
})


cbind.BeadLevelList <- function(..., deparse.level=1) {
  #function for combining BeadLevelLists.
  #Modified from limma source code.
	objects <- list(...)
	nobjects <- length(objects)
	out <- objects[[1]]
	if(nobjects > 1)
	for (i in 2:nobjects) {
		out$R <- cbind(out$R,objects[[i]]$R)
		out$G <- cbind(out$G,objects[[i]]$G)
		out$Rb <- cbind(out$Rb,objects[[i]]$Rb)
		out$Gb <- cbind(out$Gb,objects[[i]]$Gb)
		out$weights <- cbind(out$weights,objects[[i]]$weights)
		out$targets <- cbind(out$targets,objects[[i]]$targets)
                out$x <- cbind(out$x, objects[[i]]$x)
                out$y <- cbind(out$y, objects[[i]]$y)
                out$ProbeID <- cbind(out$ProbeID, objects[[i]]$ProbeID)
                out$SAMPLE <-cbind(out$SAMPLE, objects[[i]]$SAMPLE)
	}
	out
}

cbind.BeadSummaryList <- function(..., deparse.level=1) {
  #function for combining BeadSummaryLists.
	objects <- list(...)
	nobjects <- length(objects)
	out <- objects[[1]]
	if(nobjects > 1)
	for (i in 2:nobjects) {
		out$R <- cbind(out$R,objects[[i]]$R)
		out$G <- cbind(out$G,objects[[i]]$G)
		out$Rb <- cbind(out$Rb,objects[[i]]$Rb)
		out$Gb <- cbind(out$Gb,objects[[i]]$Gb)
		out$BeadStDev <- cbind(out$BeadStDev,objects[[i]]$BeadStDev)
		out$Nobeads <- cbind(out$Nobeads,objects[[i]]$Nobeads)
                out$nooutliers <- cbind(out$nooutliers, objects[[i]]$nooutliers)
 #               out$ProbeID <- cbind(out$ProbeID, objects[[i]]$ProbeID)
                out$Detection <- cbind(out$Detection, objects[[i]]$Detection)
                out$SAMPLE <- append(out$SAMPLE, objects[[i]]$SAMPLE)
                out$SAM <- append(out$SAM, objects[[i]]$SAM)
                if(!is.null(out$other)){
                  for(j in length(out$other))
                    out$other[[j]] <- cbind(out$other[[j]], objects[[i]]$other[[j]])
                }
	}
	out
}

rbind.BeadSummaryList <- function(..., deparse.level=1) {
  #function for combining BeadSummaryLists.
	objects <- list(...)
	nobjects <- length(objects)
	out <- objects[[1]]
	if(nobjects > 1)
	for (i in 2:nobjects) {
		out$R <- rbind(out$R,objects[[i]]$R)
		out$G <- rbind(out$G,objects[[i]]$G)
		out$Rb <- rbind(out$Rb,objects[[i]]$Rb)
		out$Gb <- rbind(out$Gb,objects[[i]]$Gb)
		out$BeadStDev <- rbind(out$BeadStDev,objects[[i]]$BeadStDev)
		out$Nobeads <- rbind(out$Nobeads,objects[[i]]$Nobeads)
                out$nooutliers <- rbind(out$nooutliers, objects[[i]]$nooutliers)
                out$ProbeID <- append(out$ProbeID, objects[[i]]$ProbeID)
                out$Detection <- rbind(out$Detection, objects[[i]]$Detection)
                
                
                if(!is.null(out$other)){
                  for(j in length(out$other))
                    out$other[[j]] <- rbind(out$other[[j]], objects[[i]]$other[[j]])
                }
	}
	out
}


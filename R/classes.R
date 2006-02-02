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
                  object$beadstdev <- object$beadstdev[,j, drop=FALSE]
                  object$nobeads <- object$nobeads[,j, drop=FALSE]
                  object$probeID <- object$probeID[,j, drop=FALSE]
                  object$nooutliers <- object$nooutliers[,j, drop=FALSE]
                  if(!is.null(oc)) for(k in oc) object$other[[k]] <- object$other[[k]][,j,drop=FALSE]
                }
	else
		if (missing(j)) {
                  object$R <- object$R[i,, drop=FALSE]
                  object$G <- object$G[i,, drop=FALSE]
                  object$Rb <- object$Rb[i,, drop=FALSE]
                  object$Gb <- object$Gb[i,, drop=FALSE]
                  object$beadstdev <- object$beadstdev[i,, drop=FALSE]
                  object$nobeads <- object$nobeads[i,, drop=FALSE]
                  object$probeID <- object$probeID[i,, drop=FALSE]
                  object$nooutliers <- object$nooutliers[i,, drop=FALSE]
                  if(!is.null(oc)) for(k in oc) object$other[[k]] <- object$other[[k]][i,,drop=FALSE]
		} else {
                  object$R <- object$R[i,j, drop=FALSE]
                  object$G <- object$G[i,j, drop=FALSE]
                  object$Rb <- object$Rb[i,j, drop=FALSE]
                  object$Gb <- object$Gb[i,j, drop=FALSE]
                  object$beadstdev <- object$beadstdev[i,j, drop=FALSE]
                  object$nobeads <- object$nobeads[i,j, drop=FALSE]
                  object$probeID <- object$probeID[i,j, drop=FALSE]
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
                  object$probeID <- object$probeID[,j, drop=FALSE]
                  object$targets <- object$targets[j,,drop=FALSE]
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
                  object$probeID <- object$probeID[i,, drop=FALSE]
                  if(!is.null(oc)) for(k in oc) object$other[[k]] <- object$other[[k]][i,,drop=FALSE]
   		} else {
                  object$R <- object$R[i,j, drop=FALSE]
                  object$G <- object$G[i,j, drop=FALSE]
                  object$Rb <- object$Rb[i,j, drop=FALSE]
                  object$Gb <- object$Gb[i,j, drop=FALSE]
                  object$x <- object$x[i,j, drop=FALSE]
                  object$y <- object$y[i,j, drop=FALSE]
                  object$probeID <- object$probeID[i,j, drop=FALSE]
                  object$targets <- object$targets[j,,drop=FALSE]
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
                out$probeID <- cbind(out$probeID, objects[[i]]$probeID)
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
		out$beadstdev <- cbind(out$beadstdev,objects[[i]]$beadstdev)
		out$nobeads <- cbind(out$nobeads,objects[[i]]$nobeads)
                out$nooutliers <- cbind(out$nooutliers, objects[[i]]$nooutliers)
                out$probeID <- cbind(out$probeID, objects[[i]]$probeID)
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


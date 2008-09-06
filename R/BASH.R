## BASH functions, and other image manipulation functions:

## *** GENERATE E ***
##this function returns an E value for each bead - it finds the set of all beads on the array with the same probeID, takes the median of their intensities, and subtracts this off.
generateE <- function(BLData, array, neighbours = NULL, log = TRUE, method = "median", what = "residG", bgfilter = "none", invasions = 20)
{
	an <- arrayNames(BLData)
	method <- match.arg(method, choices = c("mean", "median"))
	bgfilter <- match.arg(bgfilter, choices = c("none", "median", "mean", "MAD", "medianMAD"))

	data <- getArrayData(BLData, array=array, log=log, trim = 0, method = method, what= what)

	#separately calculate resids for 0 probes (discarded in the createBeadSummary step within getArrayData)
	if(what %in% c("residG", "residR", "residM"))
	{
		sel <- which(BLData[[an[array]]]$ProbeID == 0)
		if(length(sel) > 0)
		{
			residwhat <- switch(EXPR = what, residG = "G", residR = "R", residM = "M")
			vals = getArrayData(BLData, array = array, log=log, what = residwhat)
			if(method == "median")
			{
				med = median(vals[which(!is.na(vals[sel]))])
				new <- c(vals[sel] - med, data[-sel])
			}
			else if(method == "mean")
			{
				valmean = mean(vals[which(!is.na(vals[sel]))])
				new <- c(vals[sel] - valmean, data[-sel])
			}
		}
		else {new <- data}
	}
	else {new <- data}

	#deal with negative or NA vals
	minG <- min(new[which(!is.na(new))])
	new <- ifelse(is.na(new), minG, new)
	E <- new

	##perform a filter afterwards if appropriate
	if(bgfilter != "none")
	{
		##generate missing data
		if(is.null(neighbours))
		{
			cat("Neighbours not specified - Generating Neighbours, using defaults...\n")
			neighbours = generateNeighbours(BLData, array)
		}
		E <- BGFilter(E = E, neighbours = neighbours, method = bgfilter, invasions = invasions)
		
	}
	E
}

## *** GENERATE NEIGHBOURS ***
generateNeighbours <- function(BLData, array, window = 30, margin = 10, thresh = 2.2)
{
	an <- arrayNames(BLData)

	data <- NULL
	data <- BLData[[an[array]]]

	wcols <- ceiling(max(data$GrnX)/(2*window))
	wrows <- ceiling(max(data$GrnY)/(2*window))

	neighbours = matrix(0, length(data$GrnX), 6)

	Coutput <- .C("Neighbours", as.double(data$GrnX), as.double(data$GrnY), as.integer(length(data$GrnX)), as.integer(t(neighbours)), as.double(thresh), as.double(margin), as.double(window), as.integer(wcols), as.integer(wrows))

	neighbours <- matrix(Coutput[[4]], ncol = 6, byrow = TRUE)
}

# *** BACKGROUND FILTERS ***
# find and subtract/scale by local background, the slow change of which is causing problems on the beadchips
BGFilter <- function(E = NULL, neighbours, invasions = 20, method = "median")
{
	##quick checks
	if(is.null(E)) {stop("No error image (E). Please run generateE to obtain one.")}
	method <- match.arg(method, choices = c("median", "mean", "MAD", "medianMAD"))
	method <- switch(method, median = 1, mean = 2, MAD = 3, medianMAD = 4)

	##C code here to apply the median filter
	Etilde <- rep(0, length(E))
	output <- .C("BGFilter", as.double(E), as.double(Etilde), as.integer(t(neighbours)), as.integer(nrow(neighbours)), as.integer(invasions), as.integer(method))
	Etilde <- output[[2]]

	##Variance: Singleton causes NAs (singletons are insignificant so no harm done)
	if(method > 2) {Etilde <- ifelse(!is.na(Etilde),Etilde,1)}

	##return 
	Etilde
}

##weighted
BGFilterWeighted <- function(E = NULL, neighbours, invasions = 20, weights = NULL)
{
	if(is.null(weights)) {weights <- rep(1, nrow(neighbours))}
	if(is.null(E)) {stop("No error image (E). Please run generateE to obtain one.")}

	##C code here to apply the filter
	Etilde <- rep(0, length(E))
	output <- .C("BGFilterWeighted", as.double(E), as.double(Etilde), as.integer(t(neighbours)), as.double(weights), as.integer(nrow(neighbours)), as.integer(invasions))
	Etilde <- output[[2]]

	##Variance: Singleton causes NAs (singletons are insignificant so no harm done)
	if(method > 2) {Etilde <- ifelse(!is.na(Etilde),Etilde,1)}

	Etilde
}

## *** CHOOSE CLUSTERS ***
## floods, returning significant clusters

chooseClusters <- function(IDs, neighbours, cutoff = 8)
{
	##first, we blank all links to non-outliers.
	neighbours2 <- matrix((neighbours %in% IDs),ncol = 6)*neighbours

	##now we can flood fill the outliers, obtaining a clusterID for each point, and a cluster size for each clusterID
	clusterID <- rep(0, nrow(neighbours))
	clustersize <- rep(0, nrow(neighbours))
	nIDs <- length(IDs)
	Coutput <- .C("FloodFill", as.integer(t(neighbours2)), as.integer(IDs), as.integer(nIDs), as.integer(clusterID), as.integer(clustersize))

	##only return beads in large enough clusters
	size <- Coutput[[5]]
	size <- size[which(size>0)]
	sig <- which(size > cutoff)
	which(Coutput[[4]] %in% sig)
}

## *** CLOSING ALGORITHM ***
## given neighbours matrix and IDs of nodes, returns the closure of these nodes

closeImage <- function(IDs, neighbours, cinvasions = 10)
{
	IDs <- sort(IDs)
	nIDs <- length(IDs)
	nbeads <- nrow(neighbours)
	IDs <- c(IDs, rep(0, nbeads - nIDs))
	output <- .C("Close", as.integer(IDs), as.integer(nIDs), as.integer(t(neighbours)), as.integer(nbeads), as.integer(cinvasions))
	sort(output[[1]][which(output[[1]] > 0)])
}

## *** DENSE REGIONS ***
##given an outlier image, return the (unclosed) locus of pixels (i.e. diffuse defects)
denseRegions <- function (IDs, neighbours, ignore = NULL, sig = 0.0001, invasions = 10)
{
	##remove ignored beads (usually compacts)
	nignore <- length(ignore)
	if(nignore > 0)
	{
		IDs <- IDs[which(!(IDs %in% ignore))]
		
		##viciously sever all links in or out of ignored beads
		keep <- 1:nrow(neighbours) %in% ignore
		keep <- !keep
		neighbours <- neighbours*keep
		neighbours <- neighbours*(!(neighbours %in% ignore))
	}

	##sort IDs
	IDs <- sort(IDs)
	nIDs <- length(IDs)
	nbeads <- nrow(neighbours)

	##call C function
	Coutput <- .C("DiffuseDefects", as.integer(t(neighbours)), as.integer(IDs), as.integer(nbeads), as.integer(nIDs), as.integer(nignore), as.integer(invasions), output = as.double(rep(0, nbeads)), as.double(sig))

	return(which(Coutput$output == 1))
}

## *** VIEW BEADS ***
##diagnostic code to take a look at the bead network
viewBeads <- function (BLData, array, x, y, xwidth = 100, ywidth = 100, neighbours = NULL, mark = NULL, markcol = "blue", markpch = 21, inten = TRUE, low = "black", high = "green", what = "G", log = TRUE, zlim = NULL, ...)
{
	xwidth = xwidth/2
	ywidth = ywidth/2
	low = col2rgb(low)
	high = col2rgb(high)

	if((xwidth > 150 | ywidth > 150) & !is.null(neighbours))
	{
		choice = menu(c("Yes", "No", "Yes, but don't plot links"), title = "Alert: neighbours specified, and xwidth or ywidth is over 300 - this may cause severe slowdown or a crash, are you sure?")
		if(choice == 2) {stop("Cancelled.")}
		if(choice == 3) {neighbours = NULL}
	}

	an <- arrayNames(BLData)
	data <- BLData[[an[array]]]
	data <- cbind(as.numeric(row.names(BLData[[an[array]]])), data)
	data$GrnY <- max(data$GrnY) - data$GrnY
	data$G <- getArrayData(BLData, what = what, array = array, log = log, ...)
	data$G <- ifelse(is.na(data$G),min(data$G[which(!is.na(data$G))]),data$G)
	colnames(data)[1] <- "ID"

	sel <- which(abs(data$GrnX-x) < xwidth & abs(data$GrnY-y) < ywidth)

	vertex <- NULL
	vertex$x <- data$GrnX[sel]
	vertex$y <- data$GrnY[sel]
	vertex$ID <- data$ID[sel]
	if(is.null(zlim)) {vertex$colour <- (data$G[sel] - min(data$G))/(max(data$G) - min(data$G))}
	else
	{
		vertex$colour <- (data$G[sel] - zlim[1])/(zlim[2] - zlim[1])
		vertex$colour <- pmin(vertex$colour,1)
		vertex$colour <- pmax(vertex$colour,0)
	}

	plot(vertex$x, vertex$y)#, pch=19, col = rgb(0,vertex$colour,0,255,maxColorValue = 255))

	if(!is.null(neighbours))
	{
		for(i in 1:length(vertex$x))
		{
			iID <- vertex$ID[i]
			for(j in 1:6)
			{
				if (neighbours[iID,j] != 0)
				{
					lines(c(vertex$x[i],data$GrnX[neighbours[iID,j]]),c(vertex$y[i],data$GrnY[neighbours[iID,j]]))	
				}
			}
		}
	}

	colours <- rgb(low[1]+(high[1]-low[1])*vertex$colour,low[2]+(high[2]-low[2])*vertex$colour,low[3]+(high[3]-low[3])*vertex$colour,255,maxColorValue = 255)

	if(inten) {points(vertex$x, vertex$y, pch=19, col = colours)}
	if(length(mark) > 0) {points(data$GrnX[mark], data$GrnY[mark], col = markcol, pch = markpch)}
}

## log fns and outlier functions (not exported)

log2.na = function (x, ...)
{
    log2(ifelse(x > 0, x, NA), ...)
}

log2cap <- function(x)
{
	temp <- min(x[x > 0])
	x <- ifelse(x > 0, x, temp)
	x <- log2(x)
}

findAllOutliersE = function(BLData, array=1, Ehat, log = FALSE, n = 3, ignore = NULL)
{
	probeList = getArrayData(BLData, array=array, what="ProbeID")
	probes = sort(unique(probeList[probeList>0]))
	inten = Ehat ##change 1 (use Ehat, not intensities)
	nasinf = is.na(inten) | !is.finite(inten)
	nasinf = nasinf | (1:length(nasinf) %in% ignore) ##change 2 (ignore)
	inten = inten[!nasinf]
	probeList = probeList[!nasinf]
	nbeads = length(inten)
	start = 0

	foo <- .C("findAllOutliers", as.double(inten), binStatus = integer(length = nbeads), as.integer(probeList), as.integer(probes), as.integer(length(probes)), as.integer(nbeads), as.integer(start), as.double(n), PACKAGE = "beadarray")

	sel = which((probeList > 0) & (foo$binStatus == 0))
	which(!nasinf)[sel]
}

findAllOutliersIgnore <- function (BLData, array = 1, log = FALSE, n = 3, what = "G", ignore = NULL)
{
    probeList = getArrayData(BLData, array = array, what = "ProbeID")
    probes = sort(unique(probeList[probeList > 0]))
    inten = getArrayData(BLData, array = array, log = log, what = what)
    nasinf = is.na(inten) | !is.finite(inten)
	nasinf = nasinf | (1:length(nasinf) %in% ignore)
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

## *** BASH FUNCTIONS ***
## Pipeline functions

BASHCompact <- function(BLData, array, neighbours = NULL, log = TRUE, maxiter = 10, cutoff = 8, cinvasions = 10, ...)
{
	start = maxiter

	##generate missing data?
	if (is.null(neighbours))
	{
		cat("Neighbours not specified - Generating Neighbours, using defaults...\n")
		neighbours <- generateNeighbours(BLData, array)
	}

	output = NULL
	o <- 42
	while(maxiter > 0 & length(o) > 0)
	{
		o <- findAllOutliersIgnore(BLData, array = array, log = log, ignore = output, ...)
		o <- chooseClusters(o, neighbours, cutoff = cutoff)
		#o <- closeImage(o, neighbours, invasions = cinvasions)
		output <- sort(c(output,o))
		output <- closeImage(output, neighbours, cinvasions = cinvasions)
		maxiter = maxiter - 1
	}
	if(maxiter == 0){warning(paste(start,"iterations performed. (Array may be seriously defective.)"))}
	output
}

BASHDiffuse <- function(BLData, array, neighbours = NULL, E = NULL, n = 3, compact = NULL, sig = 0.0001, invasions = 10, cutoff = 8, cinvasions = 10, twotail = FALSE)
{
	##generate missing data
	if (is.null(neighbours))
	{
		cat("Neighbours not specified - Generating Neighbours, using defaults...\n")
		neighbours <- generateNeighbours(BLData, array)
	}
	if (is.null(E))
	{
		cat("Error image, E, not specified - Generating E, using bgfilter = \"median\"...\n")
		E <- generateE(BLData, array, neighbours, log = TRUE, method = "median", bgfilter = "median")
	}

	o <- findAllOutliersE(BLData, array, E, log = FALSE, n = n, ignore = compact)

	if(twotail)
	{
		##find tail outliers
		o.top <- which(E[o] > 0)
		o.bottom <- which(E[o] < 0)
		##diffuse
		diffuse.top <- denseRegions(o.top, neighbours, sig = sig, invasions = invasions, ignore = compact)
		diffuse.bottom <- denseRegions(o.bottom, neighbours, sig = sig, invasions = invasions, ignore = compact)
		diffuse <- unique(c(diffuse.top, diffuse.bottom))
	}
	else
	{
		##diffuse
		diffuse <- denseRegions(o, neighbours, sig = sig, invasions = invasions, ignore = compact)
		#diffuse <- which(diffuse > 1 - sig)
	}

	##only large clusters
	diffuse <- chooseClusters(diffuse, neighbours, cutoff = cutoff)

	##close
	diffuse <- closeImage(diffuse, neighbours, cinvasions = cinvasions)

	diffuse
}

BASHExtended <- function(BLData, array, neighbours = NULL, E = NULL, E.BG = NULL)
{
	##generate missing data
	if (is.null(neighbours) && is.null(E.BG)) ##we only need neighbours if E.BG isn't specified
	{
		cat("Neighbours not specified - Generating Neighbours, using defaults...\n")
		neighbours <- generateNeighbours(BLData, array)
	}
	if (is.null(E))
	{
		cat("Error Image, E, not specified - Generating E (with bgfilter = none)...\n")
		E <- generateE(BLData, array, neighbours, log = TRUE, method = "median", bgfilter = "none")
	}
	if (is.null(E.BG))
	{
		cat("Background Error Image, E.BG, not specified - Generating E.BG, via median filter on E.\n")
		E.BG <- E - BGFilter(E = E, neighbours = neighbours, invasions = 20, method = "median")
	}

	return(var(E.BG)/var(E))
}

## *** BASH FUNCTION ***
## Entire pipeline analysis

BASH <- function(BLData, array, compact = TRUE, diffuse = TRUE, extended = TRUE, log = TRUE, cinvasions = 10, dinvasions = 15, einvasions = 20, bgcorr = "median", maxiter = 10, compcutoff = 8, compdiscard = TRUE, diffcutoff = 10, diffsig = 0.0001, diffn = 3, difftwotail = FALSE)
{
	##checks
	bgcorr = match.arg(bgcorr, c("none", "median", "medianMAD"))
	output <- NULL
	scores <- NULL
	an <- arrayNames(BLData)
	if(is.null(array)) {array = 1:length(an)}
	if(!compact & !diffuse & !extended) {stop("All analyses disabled.")}

	for(i in array)
	{
		cat("Array",i,":\n")
		cat("Generating neighbours matrix... ")
		neighbours <- generateNeighbours(BLData, i)
		cat("done.\n")

		##COMPACT
		if(compact)
		{
			cat("Compact analysis... ")
			c <- BASHCompact(BLData, array = i, neighbours = neighbours, log = log, maxiter = maxiter, cutoff = compcutoff, cinvasions = cinvasions)
			cat("done.",length(c),"beads identified.\n")
		}
		else
		{c = numeric(0)}

		##EXTENDED
		if(extended)
		{
			cat("Extended analysis... ")
			E <- generateE(BLData, i,neighbours = neighbours, log = log, method = "median", bgfilter = "none")
			E.med <- BGFilter(E = E, neighbours = neighbours, invasions = 20, method = "median")
			E.BG <- E - E.med ##this reverses the subtraction process BGFilter performs, to get E.BG = median filtered
			scores[i] <- BASHExtended(BLData, array = i, neighbours = neighbours, E = E, E.BG = E.BG)
			cat("done. Value is",scores[i],".\n")
		}


		##DIFFUSE
		if(diffuse)
		{
			cat("Diffuse analyis... ")
			
			##if we did the extended analysis, then pass the baton to save recomputation
			if(extended)
			{
				if(bgcorr == "medianMAD") {E <- BGFilter(E = E.med, neighbours = neighbours, method = "MAD")}
				if(bgcorr == "median") {E <- E.med}
				if(bgcorr == "none") {E <- generateE(BLData, i, neighbours = neighbours, log = log, method = "median", bgfilter = bgcorr, invasions = einvasions)}
			}
			else
			{
				E <- generateE(BLData, i, neighbours = neighbours, log = log, method = "median", bgfilter = bgcorr, invasions = einvasions)
			}

			if(compdiscard)
			{
				d <- BASHDiffuse(BLData, array = i, neighbours = neighbours, E = E, n = diffn, compact = c, sig = diffsig, invasions = dinvasions, cutoff = diffcutoff, cinvasions = cinvasions, twotail = difftwotail)
			}
			else
			{
				d <- BASHDiffuse(BLData, array = i, neighbours = neighbours, E = E, n = diffn, compact = NULL, sig = diffsig, invasions = dinvasions, cutoff = diffcutoff, cinvasions = cinvasions, twotail = difftwotail)
			}
			cat("done.",length(d),"beads identified.\n")
		}
		else
		{d = numeric(0)}

		output$wts[[i]] <- as.numeric(!1:nrow(BLData[[an[i]]]) %in% unique(c(c,d)))
		cat("Weights found. Total no of defective beads:",length(which(!output$wts[[i]])),"\n")

	}
	output$ext <- scores
	output$call <- match.call()
	output
}

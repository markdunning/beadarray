#Masking GUI functions

#this function is used to find if points (x,y) are in polygon px,py - it shouldn't be exported
IsInPolygon <- function(x, y, px, py)
{
	#px, py are vectors containing the x and y co-ords of the polygon

	n <- length(px)
	stopifnot(n == length(py))
	stopifnot(n > 2)
	stopifnot(length(x) == length(y))
	stopifnot(is.numeric(x))
	stopifnot(is.numeric(y))

	nodes <- rep(FALSE, length(x))

	for(i in 1:n)
	{
		j <- i%%n + 1

		##which points are on the LHS of this line, and have y values in the same range? Flip truth of each one.
		nodes <- xor(nodes, (px[i] + (px[j] - px[i])*abs((y-py[i])/(py[j]-py[i])) > x & max(py[i], py[j]) > y & min(py[i], py[j]) <= y))
	}
	return(as.integer(nodes))
}

addArrayMask <- function(BLData, array = 0, SAM = FALSE, nrow = 50, ncol = 50, high = "red", low = "yellow", zlim = c(7,15), override = FALSE)
{
	if(class(BLData) != "BeadLevelList")
	{stop(paste("Must be performed on a BeadLevelList object."))}

	if(array == 0)
	{
		stop(paste("Please specify an array."))
	}

	vertices = NULL
	an = arrayNames(BLData)
	options(locatorBell = FALSE)

	imageplot(BLData, array = array, nrow = nrow, ncol = ncol, high = high, low = low, zlim = zlim, main = an[array])
	axis(1)
	axis(2)
	temp <- NULL
	temp <- locator(512, type = "o")
	if(length(temp$x) == 0)
		{	
			temp$x <- 0
			temp$y <- 0
		}
	vertices <- temp

	cat("Working...\n")

	if(length(vertices$x) > 2)
	{

	#rescale from imageplot to (x,y) co-ords
	vertices$x <- vertices$x / ncol
	vertices$y <- vertices$y / nrow

	x <- getArrayData(BLData, array = array, what = "GrnX")
	y <- getArrayData(BLData, array = array, what = "GrnY")
	x.min = min(x)
	x.max = max(x)
	y.min = min(y)
	y.max = max(y)
	x.range = x.max - x.min
	y.range = y.max - y.min

	vertices$x <- x.min + x.range*vertices$x
	#technical quirk (i.e. choice of origin) requires flipping of polygon for use with raw BL data
	vertices$y <- y.max - y.range*vertices$y

	##find selected beads using IsInPolygon
	Masked <- IsInPolygon(x, y, vertices$x, vertices$y)
	wts <- 1 - Masked ##0 in this vector indicates new mask to apply

	##warn here if you will eliminate a probe type
	probes <- BLData[[an[array]]]$ProbeID
	before <- sort(unique(probes))
	if(is.null(BLData[[an[array]]]$wts))
	{
		after <- probes[which(wts == 1)]
	}
	else
	{
		sel <- which(wts*BLData[[an[array]]]$wts == 1)
		after <- BLData[[an[array]]]$ProbeID[sel]
	}
	after <- sort(unique(after))
	
	if(length(before) != length(after))
	{
		cat("WARNING: This mask would eliminate", length(before) - length(after), "probe type(s).")
		if(length(before) - length(after) <= 500) {cat(" Displaying locations on the plot as crosses.")}
		cat("\n")
	}

	#display (unless there are too many)
	o <- findAllOutliers(BLData, array = array)
	plotBeadLocations(BLData, array = array, BeadIDs = o, main = an[array], SAM = SAM, pch = ".")

	if(!is.null(BLData[[an[array]]]$wts))
	{
		sel <- which(BLData[[an[array]]]$wts == 0)
		if(override | length(sel) <= 200000)
		{
			points(x[sel], y.max - y[sel], pch = ".", col = "Grey")
		}
		else {cat("WARNING: Current mask not plotted, since there are over 200 000 beads in the mask. (You can override this by setting override = TRUE.)\n")}
	}
	polygon(vertices$x, y.max - vertices$y, density = c(10,20), col = "Red")
	if(length(before) - length(after) <= 500 & length(before) - length(after) > 0 )
	{
		##display eliminated beads
		eliminatedIDs <- before[which(!(before %in% after))]
		sel <- which(probes %in% eliminatedIDs)
		points(x[sel], y.max - y[sel], pch = "X", col = "Blue")
	}

	##"Is this ok?" prompt
	choice <- menu(c("Accept", "Accept and define another region", "Reject", "Reject and redefine"), title = "Accept new region for masking?")

	if(choice == 1 | choice == 2)
	{
		##now makes changes permanent
		##either add column or alter column
		if(is.null(BLData[[an[array]]]$wts))
		{
			#wts <- as.data.frame(wts)
			assign(an[array], cbind(BLData[[an[array]]], wts) , envir = BLData@beadData)
		}
		else
		{
			wts <- pmin(wts,BLData[[an[array]]]$wts)
			temp <- BLData[[an[array]]]
			temp$wts <- wts
			assign(an[array], temp , envir = BLData@beadData)
		}
	}

	if(choice == 2 | choice == 4)
	{
		addArrayMask(BLData, array, SAM, nrow, ncol, high, low, zlim)
	}

	}

	else{cat("Not a polygon.")}
}

removeArrayMask <- function(BLData, array = 0, SAM = FALSE, nrow = 50, ncol = 50, high = "red", low = "yellow", zlim = c(7,15), override = FALSE)
{
	if(class(BLData) != "BeadLevelList")
	{stop(paste("Must be performed on a BeadLevelList object."))}

	if(array == 0)
	{
		stop(paste("Please specify an array."))
	}

	vertices = NULL
	an = arrayNames(BLData)
	options(locatorBell = FALSE)

	##Return null if no weights
	if(is.null(BLData[[an[array]]]$wts))
	{stop(paste("No weights on BeadLevelList object. Please define weights using e.g. addArrayMask first."))}

	imageplot(BLData, array = array, nrow = nrow, ncol = ncol, high = high, low = low, zlim = zlim, main = an[array])
	axis(1)
	axis(2)
	temp <- NULL
	temp <- locator(512, type = "o")
	if(length(temp$x) == 0)
	{	
		temp$x <- 0
		temp$y <- 0
	}
	vertices <- temp

	cat("Working...\n")

	if(length(vertices$x) > 2)
	{
		#rescale from imageplot to (x,y) co-ords
		vertices$x <- vertices$x / ncol
		vertices$y <- vertices$y / nrow
	
		x <- getArrayData(BLData, array = array, what = "GrnX")
		y <- getArrayData(BLData, array = array, what = "GrnY")
		x.min = min(x)
		x.max = max(x)
		y.min = min(y)
		y.max = max(y)
		x.range = x.max - x.min
		y.range = y.max - y.min

		vertices$x <- x.min + x.range*vertices$x
		#technical quirk (i.e. choice of origin) requires flipping of polygon for use with raw BL data
		vertices$y <- y.max - y.range*vertices$y

		##find selected beads using IsInPolygon
		Masked <- IsInPolygon(x, y, vertices$x, vertices$y)
		wts <- Masked ## 1 in this vector indicates an area to unmask

		##warn here if you will still eliminate a probe type
		probes <- BLData[[an[array]]]$ProbeID
		before <- sort(unique(probes))

		sel <- which(BLData[[an[array]]]$wts + wts != 0) ##which(pmax(BLData[[an[array]]]$wts,wts) == 1)
		after <- probes[sel]
		after <- sort(unique(after))
		
		if(length(before) != length(after))
		{
			cat("WARNING: This mask would still eliminate", length(before) - length(after), "probe type(s).")
			if(length(before) - length(after) <= 500) {cat(" Displaying locations on the plot as crosses.")}
			cat("\n")
		}

		#display
		o <- findAllOutliers(BLData, array = array)
		plotBeadLocations(BLData, array = array, BeadIDs = o, main = an[array], SAM = SAM, pch = ".")

		if(!is.null(BLData[[an[array]]]$wts))
		{
			sel <- which(BLData[[an[array]]]$wts == 0)
			if(override | length(sel) <= 200000)
			{
				points(x[sel], y.max - y[sel], pch = ".", col = "Grey")
			}
			else {cat("WARNING: Current mask not plotted, since there are over 200 000 beads in the mask. (You can override this option by setting override = TRUE.)\n")}
		}
		polygon(vertices$x, y.max - vertices$y, density = c(10,20), col = "Red")

		if(length(before) - length(after) <= 500 & length(before) - length(after) > 0)
		{
			##display eliminated beads
			eliminatedIDs <- before[which(!(before %in% after))]
			sel <- which(BLData[[an[array]]]$ProbeID %in% eliminatedIDs)
			points(x[sel], y.max - y[sel], pch = "X", col = "Blue")
		}

		##"Is this ok?" prompt
		choice <- menu(c("Accept", "Accept and define another region", "Reject", "Reject and redefine"), title = "Accept new region for unmasking?")

		if(choice == 1 | choice == 2)
		{		
			##Make changes permanent
			wts <- pmax(BLData[[an[array]]]$wts,wts)
			temp <- BLData[[an[array]]]
			temp$wts <- wts
			assign(an[array], temp , envir = BLData@beadData)
			if(length(BLData[[an[array]]]$wts[BLData[[an[array]]]$wts != 1]) == 0) 
			{
				cat("Mask entirely removed. Clearing weights.\n")
				clearArrayMask(BLData, array)
			}

		}

		if(choice == 2 | choice == 4)
		{
			removeArrayMask(BLData, array, SAM, nrow, ncol, high, low, zlim)
		}

	}
	else{cat("Not a polygon.")
	}
}

showArrayMask <- function(BLData,array = 0, SAM = FALSE, elim = TRUE, override = FALSE)
{
	if(class(BLData) != "BeadLevelList")
	{stop(paste("Must be performed on a BeadLevelList object."))}

	if(array == 0)
	{
		stop(paste("Please specify an array."))
	}

	an = arrayNames(BLData)

	if(!is.null(BLData[[an[array]]]$wts))
	{
		sel <- which(BLData[[an[array]]]$wts==0)

		##get out now if there are too many masked beads
		if(!override & length(sel) > 200000)
		{
			stop(paste("There are over 200 000 beads in the mask, thus plotting the mask may cause R to freeze. (You can override this error by setting override = TRUE.)\nNumber of masked beads:",length(sel),"\n"))
		}

		o <- findAllOutliers(BLData, array = array)
		plotBeadLocations(BLData, array = array, BeadIDs = o, main = an[array], SAM = SAM, pch = ".")

		y.max = max(getArrayData(BLData, array = array, what = "GrnY"))

		x.cds <- BLData[[an[array]]]$GrnX[sel]
		y.cds <- BLData[[an[array]]]$GrnY[sel]
		points(x.cds, y.max - y.cds, pch = ".", col = "Red")

		if(elim)
		{
			eliminated <- listEliminatedProbes(BLData, 1)
			sel <- which(BLData[[an[array]]]$ProbeID %in% eliminated)
			x.cds <- BLData[[an[array]]]$GrnX[sel]
			y.cds <- BLData[[an[array]]]$GrnY[sel]
			points(x.cds, y.max - y.cds, pch = "X", col = "Blue")			
		}
	}
}

listEliminatedProbes <- function(BLData, array = 0)
{
	if(class(BLData) != "BeadLevelList")
	{stop(paste("Must be performed on a BeadLevelList object."))}

	if(array == 0)
	{
		stop(paste("Please specify an array."))
	}

	an = arrayNames(BLData)

	if(is.null(BLData[[an[array]]]$wts))
	{
		return(vector("numeric",0))
	}
	else
	{
		before <- sort(unique(BLData[[an[array]]]$ProbeID))
		sel <- which(BLData[[an[array]]]$wts == 1)
		after <- BLData[[an[array]]]$ProbeID[sel]
		after <- sort(unique(after))
		before[which(!(before %in% after))]
	}
}

clearArrayMask <- function(BLData,array = 0)
{
	if(class(BLData) != "BeadLevelList")
	{stop(paste("Must be performed on a BeadLevelList object."))}

	if(array == 0)
	{
		stop(paste("Please specify an array."))
	}

	an = arrayNames(BLData)

	if(!is.null(BLData[[an[array]]]$wts))
	{
		BLData <- setWeights(BLData, array = array, wts = 1)
		cat("Mask removed.\n")
	}
	else
	{
		cat("No mask to remove.\n")
	}
}

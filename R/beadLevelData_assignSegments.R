assignSegments = function(BLData, array=1, nSegments = 8, useLocs = F){

	ys = getBeadData(BLData, what="GrnY", array=array)

	ys = ys - min(ys)

	segEnds = seq(from=0, to = max(ys), by = max(ys)/9)

	cut(ys, segEnds)

}
